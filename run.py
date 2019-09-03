#!/usr/bin/env python3

import datetime
import json
import logging
import multiprocessing
import os
import pprint
import time

import dateutil.parser
import flywheel
import flywheel.rest
from google.api_core.exceptions import NotFound
from google.auth.transport.requests import Request
from google.cloud import bigquery, storage
from google.oauth2 import service_account
from google.oauth2.credentials import Credentials

logging.basicConfig(
    level=logging.ERROR,
)
log = logging.getLogger('genomics_deepvariant')


def main(context):
    config = context.config
    log.setLevel(getattr(logging, config['log_level']))

    # setup service account
    # TODO: add support for user accounts
    service_account_path = context.get_input_path('service_account')
    with open(service_account_path, 'r', encoding='utf-8') as json_fi:
        credentials_info = json.load(json_fi)
    credentials = service_account.Credentials.from_service_account_info(
        credentials_info,
        scopes=[
            'https://www.googleapis.com/auth/cloud-platform',
            'https://www.googleapis.com/auth/genomics'
        ]
    )

    gcp_project = config.get('project_id')

    # Initialize Google API clients
    storage_client = storage.Client(project=gcp_project, credentials=credentials)
    genomics_client = GenomicsClient(credentials)
    bq_client = bigquery.Client(project=gcp_project, credentials=credentials)

    pipeline_regions = config.get('regions', '').split(' ')
    pipeline_data_disk_size = config.get('base_pipeline_disk_size')

    dv_runner_zones = config.get('zones', '')
    dv_runner_image_version = config.get('image_version')
    dv_runner_docker_image = 'gcr.io/deepvariant-docker/deepvariant:0.7.2'.format(dv_runner_image_version)
    dv_runner_model = 'gs://deepvariant/models/DeepVariant/0.7.2/DeepVariant-inception_v3-0.7.2+data-wgs_standard'

    working_bucket_name = '{}-deepvariant-{}'.format(gcp_project, str(int(time.time())))
    # create the bucket
    storage_client.create_bucket(working_bucket_name)

    working_bucket_path = 'gs://{}'.format(working_bucket_name)
    output_filename = 'output.vcf'
    output_file_storage_path = '{}/{}'.format(working_bucket_path, output_filename)
    staging_storage_path = '{}/staging'.format(working_bucket_path)

    # upload inputs to the working bucket
    uploaded_inputs = upload_input_files(context, credentials, working_bucket_name)

    # prepare pipeline main resurce
    resources = {
        'virtualMachine': {
            'machineType': config.get('machine_type'),
            'serviceAccount': {
                'scopes': [
                    'https://www.googleapis.com/auth/compute',  # start worker machines
                    'https://www.googleapis.com/auth/devstorage.full_control',  # accessing and storing input/output files
                    'https://www.googleapis.com/auth/bigquery'  # export vsf to BigQuery
                ]
            },
            'bootDiskSizeGb': 15,
            'disks': [{
                'name': 'data',
                'sizeGb': pipeline_data_disk_size
            }]
        },
        'projectId': gcp_project,
        'regions': pipeline_regions
    }

    # start preparing actions
    actions = []

    # create the destination dataset
    bq_dataset_name = config.get('dest_bq_dataset')
    try:
        bq_client.get_dataset(bq_dataset_name)
    except NotFound:
        log.info('Creating destination BigQuery dataset... %s', bq_dataset_name)
        bq_client.create_dataset(bq_dataset_name)

    # create index of of input bam file if it doesn't exists
    if not uploaded_inputs.get('bai'):
        # donwload the file to the data disk
        actions.append({
            'imageUri': 'google/cloud-sdk',
            'commands': [
                'sh',
                '-c',
                'gsutil cp {} /data/{}'.format(
                    uploaded_inputs['bam']['path'],
                    uploaded_inputs['bam']['file_name']
                )
            ],
            'mounts': [
                {
                    'disk': 'data',
                    'path': '/data'
                }
            ]
        })
        actions.append({
            'imageUri': 'gcr.io/genomics-tools/samtools',
            'commands': [
                'index',
                '/data/' + uploaded_inputs['bam']['file_name'],
                '/data/' + uploaded_inputs['bam']['file_name'] + '.bai'
            ],
            'mounts': [
                {
                    'disk': 'data',
                    'path': '/data'
                }
            ]
        })

        uploaded_inputs['bai'] = {
            'file_name': uploaded_inputs['bam']['file_name'] + '.bai',
            'path': '{}/input/{}.bai'.format(
                working_bucket_path,
                uploaded_inputs['bam']['file_name']
            )
        }

    # donwload the file to the data disk
    actions.append({
        'imageUri': 'google/cloud-sdk',
        'commands': [
            'sh',
            '-c',
            'gsutil cp {} /data/{}'.format(
                uploaded_inputs['fasta']['path'], uploaded_inputs['fasta']['file_name']
            )
        ],
        'mounts': [
            {
                'disk': 'data',
                'path': '/data'
            }
        ]
    })

    if uploaded_inputs['fasta']['file_name'].endswith('.gz'):
        # create the index file
        actions.append({
            'imageUri': 'gcr.io/genomics-tools/samtools',
            'entrypoint': '/bin/sh',
            'commands': [
                '-c',
                'gzip -d /data/' + uploaded_inputs['fasta']['file_name']
            ],
            'mounts': [
                {
                    'disk': 'data',
                    'path': '/data'
                }
            ]
        })

    if not uploaded_inputs.get('fai'):
        # index if it was not provided
        actions.append({
            'imageUri': 'gcr.io/genomics-tools/samtools',
            'entrypoint': '/bin/sh',
            'commands': [
                '-c',
                'samtools faidx /data/{0}'.format(
                    uploaded_inputs['fasta']['file_name'].rstrip('.gz')
                )
            ],
            'mounts': [
                {
                    'disk': 'data',
                    'path': '/data'
                }
            ]
        })
    else:
        # download the index file
        actions.append({
            'imageUri': 'google/cloud-sdk',
            'commands': [
                'sh',
                '-c',
                'gsutil cp {} /data/{}.fai'.format(
                    uploaded_inputs['fai']['path'], uploaded_inputs['fasta']['file_name'].rstrip('.gz')
                )
            ],
            'mounts': [
                {
                    'disk': 'data',
                    'path': '/data'
                }
            ]
        })

    # create fa.gz, fa.gz.gzi, fa.gz.fai files
    actions.append({
        'imageUri': 'gcr.io/genomics-tools/samtools',
        'entrypoint': '/bin/sh',
        'commands': [
            '-c',
            'cp /data/{0}.fai /data/{0}.gz.fai; bgzip -i /data/{0}'.format(
                uploaded_inputs['fasta']['file_name'].rstrip('.gz')
            )
        ],
        'mounts': [
            {
                'disk': 'data',
                'path': '/data'
            }
        ]
    })

    uploaded_inputs['fai'] = {
        'file_name': uploaded_inputs['fasta']['file_name'].rstrip('.gz') + '.gz.fai',
        'path': '{}/input/{}.gz.fai'.format(
            working_bucket_path,
            uploaded_inputs['fasta']['file_name'].rstrip('.gz')
        )
    }

    uploaded_inputs['fai_gzi'] = {
        'file_name': uploaded_inputs['fasta']['file_name'].rstrip('.gz') + '.gz.gzi',
        'path': '{}/input/{}.gz.gzi'.format(
            working_bucket_path,
            uploaded_inputs['fasta']['file_name'].rstrip('.gz')
        )
    }

    # upload the files back to cloud storage
    actions.append({
        'imageUri': 'google/cloud-sdk',
        'entrypoint': '/bin/sh',
        'commands': [
            '-c',
            'gsutil -m cp /data/* {}/input'.format(working_bucket_path)
        ],
        'mounts': [
            {
                'disk': 'data',
                'path': '/data'
            }
        ]
    })

    deepvariant_runner_options = {
        'project': gcp_project,
        'zones': dv_runner_zones,
        'docker_image': dv_runner_docker_image,
        'outfile': output_file_storage_path,
        'staging': staging_storage_path,
        'model': dv_runner_model,
        'bam': uploaded_inputs['bam']['path'],
        'bai': uploaded_inputs['bai']['path'],
        'ref': uploaded_inputs['fasta']['path'],
        'ref_fai': uploaded_inputs['fai']['path'],
        'ref_gzi': uploaded_inputs['fai_gzi']['path'],
        'regions': config.get('chr_regions'),
        'shards': config.get('shards'),
        'make_examples_workers': config.get('make_examples_workers'),
        'make_examples_cores_per_worker': config.get('make_examples_cores_per_worker'),
        'make_examples_ram_per_worker_gb': config.get('make_examples_ram_per_worker_gb'),
        'make_examples_disk_per_worker_gb': config.get('make_examples_disk_per_worker_gb'),
        'call_variants_workers': config.get('call_variants_workers'),
        'call_variants_cores_per_worker': config.get('call_variants_cores_per_worker'),
        'call_variants_ram_per_worker_gb': config.get('call_variants_ram_per_worker_gb'),
        'call_variants_disk_per_worker_gb': config.get('call_variants_disk_per_worker_gb'),
        'postprocess_variants_cores': config.get('postprocess_variants_cores'),
        'postprocess_variants_ram_gb': config.get('postprocess_variants_ram_gb'),
        'postprocess_variants_disk_gb': config.get('postprocess_variants_disk_gb'),
        'preemptible': '--preemptible' if config.get('preemptible') else '',
        'max_preemptible_tries': config.get('max_preemptible_tries'),
    }

    # add the deepvariant runner action
    actions.append({
        'imageUri': 'gcr.io/cloud-genomics-pipelines/gcp-deepvariant-runner',
        'commands': [
            'bash',
            '-c',
            '/opt/deepvariant_runner/bin/gcp_deepvariant_runner '
            '--project {project} '
            '--zones {zones} '
            '--docker_image {docker_image} '
            '--outfile {outfile} '
            '--staging {staging} '
            '--model {model} '
            '--bam {bam} '
            '--bai {bai} '
            '--ref {ref} '
            '--ref_fai {ref_fai} '
            '--ref_gzi {ref_gzi} '
            '--shards {shards} '
            '--make_examples_workers {make_examples_workers} '
            '--make_examples_cores_per_worker {make_examples_cores_per_worker} '
            '--make_examples_ram_per_worker_gb {make_examples_ram_per_worker_gb} '
            '--make_examples_disk_per_worker_gb {make_examples_disk_per_worker_gb} '
            '--call_variants_workers {call_variants_workers} '
            '--call_variants_cores_per_worker {call_variants_cores_per_worker} '
            '--call_variants_ram_per_worker_gb {call_variants_ram_per_worker_gb} '
            '--call_variants_disk_per_worker_gb {call_variants_disk_per_worker_gb} '
            '--postprocess_variants_cores {postprocess_variants_cores} '
            '--postprocess_variants_ram_gb {postprocess_variants_ram_gb} '
            '--postprocess_variants_disk_gb {postprocess_variants_disk_gb} '
            '{preemptible} '
            '--max_preemptible_tries {max_preemptible_tries} '
            '--regions "{regions}" '
            '--gcsfuse'.format(
                **deepvariant_runner_options
            )
        ]
    })

    output_table = '{}:{}.vcf_export_{}'.format(gcp_project, bq_dataset_name, datetime.datetime.now().strftime('%Y%m%d_%H%M%S'))

    # add vcf to bigquery export action
    actions.append({
        'imageUri': 'gcr.io/gcp-variant-transforms/gcp-variant-transforms',
        'entrypoint': '/opt/gcp_variant_transforms/bin/vcf_to_bq',
        'commands': [
            '--input_pattern',
            output_file_storage_path,
            '--output_table',
            output_table,
            '--temp_location',
            staging_storage_path,
            '--job_name',
            'vcf-to-bq',
            '--runner',
            'DataflowRunner',
            '--project',
            gcp_project
        ]
    })

    try:
        op = genomics_client.run(
            actions,
            resources
        )
        genomics_client.wait_operation(op)
    except GCPError as e:
        log.exception(e)
        exit(1)
    else:
        log.info('Deepvariant finished')
        download_output_files(context, credentials, uploaded_inputs, output_file_storage_path)
    finally:
        log.info('Removing storage bucket: %s', working_bucket_name)
        bucket = storage_client.get_bucket(working_bucket_name)
        bucket.delete(force=True)



def upload_input_files(context, credentials, bucket_name):
    files_to_upload = [{'input_path': context.get_input_path(input_name), 'name': input_name} for input_name in ['bam', 'bai', 'fasta', 'fai']]
    files_to_upload = filter(lambda f: f['input_path'], files_to_upload)

    if not files_to_upload:
        return

    gcp_project = context.config.get('project_id')
    storage_client = storage.Client(project=gcp_project, credentials=credentials)
    pool = multiprocessing.Pool()
    uploaded_files = pool.map(upload_file, [(credentials, gcp_project, 'gs://{}/input'.format(bucket_name), f) for f in files_to_upload])
    return {f['name']: f for f in uploaded_files}


def upload_file(args_tuple):
    """Upload a single file to object storage."""
    credentials, project, gcs_prefix, input_to_upload = args_tuple
    file_path = input_to_upload['input_path']
    log.info('Uploading file %s to %s ...', file_path, gcs_prefix)
    bucket_name, _, prefix = gcs_prefix.replace('gs://', '').partition('/')
    storage_client = storage.Client(project=project, credentials=credentials)
    bucket = storage_client.get_bucket(bucket_name)
    object_name = prefix + '/' + os.path.basename(file_path).lstrip('/')
    blob = bucket.blob(object_name)
    blob.upload_from_filename(filename=file_path)
    log.info('Uploading file %s to %s ...done', file_path, gcs_prefix)
    return {
        'name': input_to_upload['name'],
        'path': 'gs://{}/{}'.format(bucket_name, object_name),
        'file_name': os.path.basename(file_path).lstrip('/')
    }


def download_output_files(context, credentials, uploaded_inputs, output_file_storage_path):
    log.info('Collecting output files ...')
    gcp_project = context.config.get('project_id')
    save_as_vcf = uploaded_inputs['bam']['file_name'].replace('.bam', '') + datetime.datetime.now().strftime('%Y%m%d_%H%M%S') + '.vcf'
    files_to_download = [(credentials, gcp_project, output_file_storage_path, save_as_vcf)]
    context.update_file_metadata(save_as_vcf, type='text')

    for input_name in ['bai', 'fai']:
        if not context.get_input_path(input_name):
            files_to_download.append((
                credentials, gcp_project, uploaded_inputs[input_name]['path'], uploaded_inputs[input_name]['file_name']
            ))
            context.update_file_metadata(uploaded_inputs[input_name]['file_name'], type='text')

    pool = multiprocessing.Pool()
    uploaded_files = pool.map(download_file, files_to_download)


def download_file(args_tuple):
    credentials, project, storage_path, save_as = args_tuple
    log.info('Downloading output file %s ...', storage_path)
    bucket_name, object_name = storage_path.replace('gs://', '').split('/', 1)
    storage_client = storage.Client(project=project, credentials=credentials)
    bucket = storage_client.get_bucket(bucket_name)
    blob = bucket.blob(object_name)
    with context.open_output(save_as, 'wb') as f:
        blob.download_to_file(f)
    log.info('Downloading output file %s ...done', storage_path)


class GenomicsClient:
    baseurl = 'https://genomics.googleapis.com'

    def __init__(self, credentials, version='v2alpha1'):
        self._credentials = credentials
        self._version = version
        self._request = Request()

    def request(self, *args, **kwargs):
        headers = kwargs.get('headers', {})
        self._credentials.before_request(self._request, *args, headers=headers)
        kwargs['headers'] = headers
        return self._request(*args, **kwargs)


    def run(self, actions, resources, environment=None, timeout=None, labels=None):
        payload = {
            'labels': labels,
            'pipeline': {
                'actions': actions,
                'resources': resources,
                'environment': environment,
                'timeout': timeout
            }
        }
        body = json.dumps(payload)
        resp = self.request('{}/{}/pipelines:run'.format(self.baseurl, self._version), 'POST', body=body)
        return json.loads(resp.data)

    def wait_operation(self, operation, sleep=10, timeout=None):
        """Wait for an operation to finish."""
        start = time.time()
        last_check = None
        while not operation['done']:
            if timeout is not None and time.time() > start + timeout:
                message = 'Wait timeout: exceeded {}s for {}'.format(timeout, operation['name'])
                raise GCPError(response=operation, message=message)

            if operation.get('error'):
                message = 'Pipeline failed'
                raise GCPError(response=operation, message=message)
            resp = self.request('{}/{}/{}'.format(self.baseurl, self._version, operation['name']), 'GET')
            operation = json.loads(resp.data)
            for event in operation['metadata']['events'][::-1]:
                timestamp = dateutil.parser.parse(event['timestamp'])
                if not last_check or timestamp > last_check:
                    log.info('\t%s\t%s', timestamp.strftime('%Y-%m-%d %H:%M:%S'), event['description'])
                    last_check = timestamp
            time.sleep(sleep)
        return operation


class GCPError(Exception):
    """GCPError class."""

    def __init__(self, response=None, message=None, status_code=None):
        """Initialize a GCPError instance."""
        if response == message is None:
            raise TypeError('response or message required')
        self.status_code = status_code or (response.status_code if response else None)
        self.message = message or self.get_message(response)
        if self.status_code is not None:
            self.message = '{}: {}'.format(self.status_code, self.message)
        super().__init__(self.message)

    @staticmethod
    def get_message(response):
        """Get error message of a respone."""
        try:
            return response.json()['error']['message']
        except Exception:  # pylint: disable=broad-except
            return response.content


if __name__ == '__main__':
    with flywheel.GearContext() as context:
        context.init_logging()
        context.log_config()
        main(context)

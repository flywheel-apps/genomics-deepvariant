FROM python:3.7.4-alpine3.10

RUN apk add --no-cache bash git

WORKDIR /flywheel/v0
COPY Pipfile Pipfile
COPY Pipfile.lock Pipfile.lock
RUN pip3 install pipenv \
    && pipenv install --deploy --system

COPY . .

CMD ["/flywheel/v0/run.py"]

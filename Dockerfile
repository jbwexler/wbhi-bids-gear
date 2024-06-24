FROM python:3.10-slim

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        ca-certificates \
        netbase \
        git \
        && \
    apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

ENV FLYWHEEL=/flywheel/v0
COPY requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir -r /tmp/requirements.txt

RUN mkdir -p ${FLYWHEEL}
COPY run.py ${FLYWHEEL}/run.py
WORKDIR ${FLYWHEEL}

ENTRYPOINT ["/flywheel/v0/run.py"]

FROM python:3.11-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    openbabel \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir poetry==2.2.1

WORKDIR /app

COPY pyproject.toml poetry.lock ./
RUN poetry config virtualenvs.create false \
    && poetry install --no-interaction --no-ansi --no-root

COPY stride /usr/local/bin/stride
RUN chmod +x /usr/local/bin/stride

COPY *.py ./

RUN mkdir -p /pdb /data /results /uniprot_cache /tmp/dbenrich26

ENTRYPOINT ["python"]
CMD ["processor.py"]

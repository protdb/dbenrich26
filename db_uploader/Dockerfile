FROM python:3.11-slim

WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY upload.py .

RUN mkdir -p /data /data/done

ENTRYPOINT ["python", "-u", "upload.py"]

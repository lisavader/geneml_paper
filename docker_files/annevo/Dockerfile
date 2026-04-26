FROM python:3.11-slim
RUN apt-get update && apt-get install -y git gcc && rm -rf /var/lib/apt/lists/*
RUN pip install setuptools
RUN git clone --branch v2.2.2 https://github.com/xjtu-omics/ANNEVO.git /annevo
RUN pip install --no-build-isolation -r /annevo/requirements.txt

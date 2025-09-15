FROM python:3.11

RUN apt-get update \
 && apt-get install -y --no-install-recommends pigz \
 && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache fire pandas numpy scipy numba matplotlib dnaio cutadapt scanpy

RUN pip install git+https://github.com/cakirb/cbutools.git@faster-firststep
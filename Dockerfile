FROM python:3.11

RUN pip install --no-cache fire pandas numpy scipy numba matplotlib dnaio cutadapt scanpy

COPY cbutools /opt/cbutools
FROM python:3.13
LABEL authors="asherpreskasteinberg"

WORKDIR /usr/src/app

COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

COPY . .
RUN pip install .

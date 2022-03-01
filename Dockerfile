# Use a dedicated base image for atat_backend
FROM ubuntu:latest

# Add basic metadata
MAINTAINER Shan Tharanga "stwm2@student.london.ac.uk"

# Update the package manager
RUN apt update -y && \
    apt install -y software-properties-common && \
    add-apt-repository ppa:deadsnakes/ppa && \
    apt update -y && \
    apt upgrade -y && \
    apt install -y python3.10 python3-pip python3.10-venv python3.10-distutils curl && \
    rm -rf /var/lib/apt/lists/*

# Then we install Poetry
ENV POETRY_HOME=/poetry
ENV PATH=$POETRY_HOME/bin:$PATH
RUN curl -sSL https://install.python-poetry.org | python3.10

WORKDIR /home/atat_backend/

# We then add the project requirement and metadata into the image
COPY poetry.toml pyproject.toml wsgi.py /home/atat_backend/

# We then install the dependancies
RUN poetry install

# Then we add the most dynamic parts of the project in to the image
COPY poetry.lock /home/atat_backend/
COPY atat_backend /home/atat_backend/atat_backend/

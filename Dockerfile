# Use a dedicated base image for viva_atat
FROM ubuntu:latest

# Add basic metadata
MAINTAINER Shan Tharanga "stwm2@student.london.ac.uk"

# Update the package manager
RUN apt update && \
    apt install -y software-properties-common && \
    add-apt-repository ppa:deadsnakes/ppa && \
    apt update && \
    apt install -y python3.9 && \
    apt install -y python3-pip && \
    apt install -y curl && \
    rm -rf /var/lib/apt/lists/*

# Then we install Poetry
ENV POETRY_HOME=/poetry
ENV PATH=$POETRY_HOME/bin:$PATH
RUN curl -sSL https://raw.githubusercontent.com/python-poetry/poetry/master/install-poetry.py | python3.9 - -p

WORKDIR /home/viva_atat/

# We then add the project requirement and metadata into the image
COPY poetry.toml pyproject.toml wsgi.py /home/viva_atat/

# We then install the dependancies
RUN poetry install

# Then we add the most dynamic parts of the project in to the image
COPY poetry.lock /home/viva_atat/
COPY viva_atat /home/viva_atat/viva_atat/

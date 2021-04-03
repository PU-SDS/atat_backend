# We will use the latest Ubuntu image
FROM ubuntu:latest

# Add basic metadata
MAINTAINER Shan Tharanga "stwm2@student.london.ac.uk"

# Then copy the repository that we checked out into the image
COPY . /home/atat

# Update the package manager
RUN apt update

# Install Python 3.8
RUN apt install software-properties-common &&\
    add-apt-repository ppa:deadsnakes/ppa &&\
    apt update &&\
    apt install python3.8
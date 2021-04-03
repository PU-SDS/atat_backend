# We will use the latest Ubuntu image
FROM ubuntu:latest

# Add basic metadata
MAINTAINER Shan Tharanga "stwm2@student.london.ac.uk"

# Then copy the repository that we checked out into the image
COPY . /home/atat

# Update the package manager
RUN apt update

# Install Python 3.8 dependancies
RUN apt install -y software-properties-common

# Add latest repository
RUN add-apt-repository ppa:deadsnakes/ppa

# Update package manager
RUN apt update

# Install Python 3.8
RUN apt install -y python3.8

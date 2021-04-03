# We will use the latest Ubuntu image
FROM ubuntu:latest

# Add basic metadata
MAINTAINER Shan Tharanga "stwm2@student.london.ac.uk"

# Then copy the repository that we checked out into the image
COPY . /home/atat

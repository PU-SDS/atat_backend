# Use a dedicated base image for viva_atat
FROM bilsab/atat-docker-base:0.1.0

# Add basic metadata
MAINTAINER Shan Tharanga "stwm2@student.london.ac.uk"

# We then add the package into the image
COPY viva_atat /home/atat/atat/

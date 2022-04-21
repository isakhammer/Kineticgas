# Importing minimal macOS container from https://hub.docker.com/r/sickcodes/docker-osx
FROM ubuntu:20.04
SHELL ["/bin/bash", "-c"]
ENV HOME_DIR /root

# ESSENTIAL
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update
RUN apt-get -y upgrade
RUN apt-get install -yq apt-utils dialog
RUN apt-get install -yq build-essential software-properties-common
RUN apt-get update

# potential dependencies
RUN apt-get -y install curl
RUN apt-get -y install gcc g++ make
RUN apt-get install -yq git cmake make

# python3.8
RUN add-apt-repository -y ppa:deadsnakes/ppa
RUN apt-get -y update
RUN apt-get install -y python3.8 python3.8-tk python3.8-dev
RUN apt-get -y install python3-pip
RUN /usr/bin/python3.8 -m pip install --upgrade pip
RUN pip3 install --upgrade pip


ENV HOME_DIR /root
RUN mkdir -p $HOME_DIR/code
ENV CODE_DIR $HOME_DIR/code
WORKDIR $CODE_DIR
# Copy all files from host to container
COPY . .

RUN pip3 install -r requirements.txt
RUN . build_dist.sh




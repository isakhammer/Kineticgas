#!/usr/bin/env bash

# Build docker image and tag "kinetic_gas:latest"
# If the build does succeed it will store the installation image as a binary
# You can check the available images on your computer using "docker image ls"
docker build . -t kinetic_gas:latest




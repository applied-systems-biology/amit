#!/bin/bash

if [[ $UID != 0 ]]; then
    echo "Please run this script with root permissions:"
    echo "sudo $0 $*"
    exit 1
fi

IMAGE_NAME="amit3"
docker build -t $IMAGE_NAME .
CONTAINER_ID=$(docker create $IMAGE_NAME)

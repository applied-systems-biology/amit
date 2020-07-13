#!/bin/bash

if [[ $UID != 0 ]]; then
    echo "Please run this script with root permissions:"
    echo "sudo $0 $*"
    exit 1
fi

IMAGE_NAME="amit3"

docker run --user 1001:1001 --volume /data/1_Projects/Tracking/:/data -it $IMAGE_NAME /bin/bash

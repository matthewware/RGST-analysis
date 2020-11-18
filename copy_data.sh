#!/bin/bash

# copy files from docker image to computer
ID=$(docker ps -q)
docker cp $ID:/home/jovyan/work/things-to-copy/ things-to-copy/


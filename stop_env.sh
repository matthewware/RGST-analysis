#!/bin/bash

# stop the analysis container
ID=$(docker ps -q --filter "ancestor=rgst_env")
docker stop $ID

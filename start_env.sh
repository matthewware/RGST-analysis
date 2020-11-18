#!/bin/bash

# start the container
docker run -P rgst_env &

sleep 3

ID=$(docker ps -q)
PORT=$(docker port $ID | awk -F '[:]+' '{ print $2 }')
URL=$(docker logs --tail 1 $ID 2>&1 | sed "s/8888/${PORT}/g")

echo 'ID: ' $ID
echo 'Port: ' $PORT
echo 'Open your broswer at: ' $URL

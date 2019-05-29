#!/usr/bin/env bash

set -x
while $1
do
 echo "Press [CTRL+C] to stop.."
 sleep 5
    echo "My second and third argument is $2 & $3"
done
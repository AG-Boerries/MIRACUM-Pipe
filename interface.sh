#!/usr/bin/env bash

# TODO: create commands with nice params
#function usage {
#    echo "usage: diz-setup_git_ssh [-ch]"
#    echo "  -c      turn on config only"
#    echo "  -h      show this help screen"
#    exit 1
#}
#
#config_only=false
#
#while getopts uh option
#do
#case "${option}"
#in
#    c) config_only='';;
#    h) usage;;
#esac
#done

set -x
while $1
do
 echo "Press [CTRL+C] to stop.."
 sleep 5
    echo "My second and third argument is $2 & $3"
done
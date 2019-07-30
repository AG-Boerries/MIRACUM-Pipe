#!/usr/bin/env bash

# TODO: create commands with nice params
function usage {
    echo "usage: diz-setup_git_ssh [-ph]"
    echo "  -p      some param"
    echo "  -h      show this help screen"
    exit 1
}

some_flag=false

while getopts ph option
do
case "${option}"
in
    p) some_flag=true;;
    h) usage;;
esac
done

echo ${some_flag}

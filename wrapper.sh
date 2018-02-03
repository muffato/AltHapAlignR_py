#!/bin/bash

# Activate the virtualenv if given a valid path
if [ ! -z "$1" ]
then
    source "$1/bin/activate"
fi
shift

exec "$@"

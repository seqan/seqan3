#!/usr/bin/env bash
GCC="/usr/bin/g++-VERSION"
DO_TIME=0

if [[ DO_TIME -eq 0 ]]; then
    $GCC "$@"
else
    FILE=$(mktemp ram_usage.XXXXXXXX)
    /usr/bin/time -v $GCC "$@" 2> $FILE
fi

#!/usr/bin/env bash

# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

# Usage: link_check.sh <SeqAn3 root directory>
# Will output the status of links in the repository.
#
# Of main interest are broken links and those with a "Link STATUS" message.
# Some URLs may not be properly matched by the regex.
#
# The general workflow is to first run the script and then check the non-working links by searching the occurrence
# within the codebase and verifying that they are indeed broken.

do_check ()
{
    RESPONSE=$(curl --http2 -Is -A 'Mozilla/5.0' $1) # HTTP2 is the default.
    if ! [[ "$RESPONSE" =~ ^HTTP.* ]]; then # If this does not work,
        RESPONSE=$(curl --http1.1 -Is -A 'Mozilla/5.0' $1) # fall back to HTTP1.1.
    fi

    HEADER=($(echo $RESPONSE | head -1)); # May look like: HTTP/2 200
    STATUS=${HEADER[1]}
    case "$STATUS" in
        200) echo "Link OK         :" $1;;
        301) echo "Link PERM MOVED :" $1;;
        302) echo "Link TEMP MOVED :" $1;;
        404) echo "Link BROKE      :" $1;;
        429) sleep 5; do_check $1;;
        *)   echo "Link STATUS" $STATUS ":" $1;;
    esac
}

if [[ $# -ne 1 ]]; then
    echo "Usage: link_check.sh <SeqAn3 root directory>"
    exit 1
fi

if [[ ! -d $1 ]]; then
    echo "The directory $1 does not exist."
    exit 1
fi

if [[ ! -f $1/include/seqan3/version.hpp ]]; then
    echo "The directory $1 does not seem to be the SeqAn3 root directory."
    echo "Cannot find $1/include/seqan3/version.hpp."
    exit 1
fi

for URL in $(grep -ohr --exclude-dir={.git,submodules} "https*://[a-zA-Z0-9./#+?=_%:-]*[a-zA-Z0-9/#+?=_%:-]" $1 | sort | uniq)
do
  do_check $URL
done

#!/usr/bin/env bash

# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0

GCC="/usr/bin/g++-VERSION"
DO_TIME=0

if [[ DO_TIME -eq 0 ]]; then
    $GCC "$@"
else
    FILE=$(mktemp ram_usage.XXXXXXXX)
    /usr/bin/time -v $GCC "$@" 2> $FILE
fi

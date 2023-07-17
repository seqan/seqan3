#!/usr/bin/env bash
# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
# -----------------------------------------------------------------------------------------------------
#
set -Eeuo pipefail

usage="\
SYNOPSIS
    update_copyright.sh <oldyear> <newyear> <file1] [<file2>] [<file3>]...

DESCRIPTION
    Updates the copyright year of files that are formatted in a certain way. And prints out the
    copyright years that it ignores.

EXAMPLES
    ./test/scripts/update_copyright.sh 2020 2021 \$(find . -not -path '*/\.*' -and -not -path '*/submodules/*' -and -not -path '*/build/*' -and -not -iname '*.patch' -type f)
        Updates all copyright entries from 2020 to 2021. Only scans non hidden directories. Does not scan build and
        submodules directory. Ignores patches in test/api_stability.
"

if [ $# -eq 0 ]; then
    echo -e "$usage"
    exit 1
fi
# New update year
oldyear=$1
year=$2
shift 2

echo "Setting copyright dates from ${oldyear} to ${year}"

for file in "$@"; do
    perl -i -pe 's/^(.*Copyright \(c\) [0-9]{4}-)'${oldyear}'(, Knut Reinert.*$)/${1}'${year}'${2}/' $file
    perl -ne 'print "'$file':$.: $_" if (/^.*Copyright.*'${oldyear}'.*$/);' $file
done

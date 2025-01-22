#!/usr/bin/env bash

# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

set -Eeuo pipefail

usage="\
SYNOPSIS
    update_copyright.sh <oldyear> <newyear> <file1> [<file2>] [<file3>]...

DESCRIPTION
    Updates the copyright year of files that are formatted in a certain way. And prints out the
    copyright years that it ignores.

EXAMPLES
    ./test/scripts/update_copyright.sh 2020 2021 \$(find . -not -path '*/\.git/**' -and -not -path '*/submodules/*' -and -not -path '*/build/*' -and -not -iname '*.patch' -type f)
        Updates all copyright entries from 2020 to 2021. Only scans non hidden directories. Does not scan build and
        submodules directory.
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
    perl -i -pe 's/^(.*SPDX-FileCopyrightText: [0-9]{4}-)'${oldyear}'(,? Knut Reinert.*$)/${1}'${year}'${2}/' $file
    perl -ne 'print "'$file':$.: $_" if (/^.*SPDX-FileCopyrightText.*'${oldyear}'.*$/);' $file
done

echo "Manually adjust the year in the argument_parser, corresponding tests, and LICENSE.md."
echo "Also provide an API Stability Patch for the argument_parser tests."

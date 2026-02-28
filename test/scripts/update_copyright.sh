#!/usr/bin/env bash

# SPDX-FileCopyrightText: 2006-2026 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2026 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

set -Eeuo pipefail

USAGE="\
SYNOPSIS
    adjust_copyright_years.sh [old_year=$(($(date +%Y) - 1))] [new_year=$(date +%Y)]

DESCRIPTION
    Updates the copyright year of files that are formatted in a certain way.

EXAMPLES
    ./test/scripts/update_copyright.sh
    ./test/scripts/update_copyright.sh 2022 # Overwrite old year
    ./test/scripts/update_copyright.sh 2024 2025 # Overwrite old and new year"

# https://www.shellcheck.net/wiki/SC2235
if [[ $# -gt 0 ]] && { [[ "$1" == "-h" ]] || [[ "$1" == "--help" ]]; }; then
    echo -e "${USAGE}"
    exit 0
fi

if [[ $# -gt 2 ]]; then
    echo -e "${USAGE}"
    exit 1
fi

OLD_YEAR="${1:-$(($(date +%Y) - 1))}"
NEW_YEAR="${2:-$(date +%Y)}"

SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
SEQAN3_ROOT="$(realpath "${SCRIPT_DIR}/../..")"

get_files()
{
    find "${SEQAN3_ROOT}" \
        -not -path "${SEQAN3_ROOT}"'/\.git/*' -and \
        -not -path "${SEQAN3_ROOT}"'/build/*' -and \
        -not -path "${SEQAN3_ROOT}"'/.vscode/*' -and \
        -not -iname "*.patch" \
        -type f \
        -print0
}

echo "Updating copyright years from ${OLD_YEAR} to ${NEW_YEAR}"

# https://www.shellcheck.net/wiki/SC2044
while IFS= read -r -d '' file; do
    perl -i -pe 's/^(.*[0-9]{4}-)'"${OLD_YEAR}"'(,? Knut Reinert.*$)/${1}'"${NEW_YEAR}"'${2}/' "${file}"
done < <(get_files)

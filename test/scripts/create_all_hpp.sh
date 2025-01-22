#!/usr/bin/env bash

# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

# Usage: create_all_hpp.sh <SeqAn3 root directory>
# Will create all.hpp in all subdirectories

# Directories (+ subdirectories) to ignore (not creating any all.hpp files in these folders)
exceptionList=(
    include/seqan3/std
    include/seqan3/contrib
    /detail
    /exposition_only
    include/seqan3/core/concept
    include/seqan3/utility/concept
)

# prints the header of a typical all.hpp file
print_header()
{
CURRENT_YEAR=$1
cat << EOF
# SPDX-FileCopyrightText: 2006-${CURRENT_YEAR} Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-${CURRENT_YEAR} Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause
EOF
}

# prints a list of "#include <file>" lines, typical for all.hpp files
print_includes()
{
    directory="$1"
    echo -ne "\n#pragma once\n\n"
    hasIncludes=0

    for path in $(find "${directory}"/* -maxdepth 0); do
        incl_path="$(realpath --relative-to include "${path}")"
        if [ "${path}" == "${directory}/all.hpp" ]; then
            continue
        elif [ -f "${path}" ]; then
            if ([[ "${path}" == *.hpp ]] || [[ "${path}" == *.h ]]) && [[ ! "${path}" == */detail.hpp ]]; then
                hasIncludes=1
               echo "#include <${incl_path}>"
            fi
        elif [[ "${path}" == */detail ]] || [[ "${path}" == */exposition_only ]]; then
            continue
        else
            # In every subdirectory (excluding detail) we generate a all.hpp file
            echo "#include <${incl_path}/all.hpp>"
            hasIncludes=1
        fi
    done
    if [[ "${hasIncludes}" -eq "0" ]]; then
        echo "#include <seqan3/core/platform.hpp>"
    fi
}

if [[ $# -ne 1 ]]; then
    echo "Usage: create_all_hpp.sh <SeqAn3 root directory>"
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
CURRENT_YEAR=$(date "+%Y")

# Iterate through all folders
for directory in $(find include/seqan3/* -type d); do
    # check if directory is in the exception list
    found=0
    for entry in "${exceptionList[@]}"; do
        if [ ! -z $(echo "${directory}" | grep $entry) ]; then
            found=1;
            break;
        fi
    done
    if [ ${found} -eq 1 ]; then
        continue
    fi

    # create all.hpp
    if [ ! -f "${directory}/all.hpp" ]; then
        echo "missing all.hpp file in ${directory}"
        (
            print_header ${CURRENT_YEAR}
            print_includes "${directory}"

        ) > "${directory}/all.hpp"
    else
        # check if license text has changed
        HEADER_LINES=$(print_header ${CURRENT_YEAR} | wc -l)
        licence_changes=$(head -n ${HEADER_LINES} "${directory}/all.hpp" | diff /dev/stdin <( print_header ${CURRENT_YEAR}) | wc -c)

        # check if include list has changed
        INCLUDE_LINES=$(print_includes "${directory}" | wc -l)
        include_changes=$(tail -n ${INCLUDE_LINES} "${directory}/all.hpp" | diff /dev/stdin <( print_includes "${directory}" ) | wc -c)

        if [ "$licence_changes" -gt 0 ] || [ "$include_changes" -gt 0 ]; then
            echo "License text or includes changed, updating ${directory}/all.hpp"
            tmp_file=$(mktemp);
            # create new license text and include list, but keep doxygen documentation
            cp "${directory}/all.hpp" "${tmp_file}"
            (
                PRAGMA_LINE=$(cat "${tmp_file}" | grep -n "#pragma once" | cut -f 1 -d ':')
                print_header ${CURRENT_YEAR}
                tail ${tmp_file} -n +$(expr 1 + ${HEADER_LINES}) | grep -v "#include" | grep -v "#pragma once"
                print_includes "${directory}"
            ) | cat -s | sed -e :a -e '/^\n*$/{$d;N;};/\n$/ba' >  "${directory}/all.hpp"
            rm $tmp_file
        fi
    fi
done

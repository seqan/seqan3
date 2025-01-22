# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

include (seqan3_path_longest_stem)

# Get a specific component of a test file which follows the seqan3 naming scheme.
# e.g. target_source_file = "range/views/take.cpp"
# component:
#  * TARGET_NAME - the target name (e.g. "take")
#  * TARGET_UNIQUE_NAME - the target name which includes the target_path (e.g. range-view-take)
#  * TEST_NAME - the test name which includes the target_path (e.g. "range/views/take")
#  * TARGET_PATH - the path to the target source (e.g. "range/view")
macro (seqan3_test_component VAR target_source_file component_name_)
    string (TOUPPER "${component_name_}" component_name)

    get_filename_component (target_relative_path "${target_source_file}" DIRECTORY)
    seqan3_path_longest_stem (target_name "${target_source_file}")

    if (component_name STREQUAL "TARGET_NAME")
        set (${VAR} "${target_name}")
    elseif (component_name MATCHES "TEST_NAME|TARGET_UNIQUE_NAME")
        if (target_relative_path)
            set (${VAR} "${target_relative_path}/${target_name}")
        else ()
            set (${VAR} "${target_name}")
        endif ()

        if (component_name STREQUAL "TARGET_UNIQUE_NAME")
            string (REPLACE "/" "-" ${VAR} "${${VAR}}")
        endif ()
    elseif (component_name STREQUAL "TARGET_PATH")
        set (${VAR} "${target_relative_path}")
    endif ()

    unset (target_name)
    unset (target_relative_path)
endmacro ()

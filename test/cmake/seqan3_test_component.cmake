# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
# -----------------------------------------------------------------------------------------------------

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
    get_filename_component (target_name "${target_source_file}" NAME_WE)

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

# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
# -----------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.10)

set (seqan3_test_include_targets
     ""
     CACHE STRING "" FORCE)

function (collect_include_target include_target)
    set (seqan3_test_include_targets
         "${seqan3_test_include_targets};${include_target}"
         CACHE STRING "" FORCE)
endfunction ()

function (list_missing_unit_tests)
    list (SORT seqan3_test_include_targets)
    foreach (include_target ${seqan3_test_include_targets})
        if (NOT TARGET ${include_target})
            continue ()
        endif ()

        get_target_property (dependencies ${include_target} MANUALLY_ADDED_DEPENDENCIES)

        # if include_target has dependencies that means a test must be associated with it
        # => so this include_target has no missing tests
        if (dependencies)
            continue ()
        endif ()

        string (REPLACE "-" "/" header ${include_target})
        string (REPLACE "_test" ".hpp" header ${header})

        get_filename_component (header_filename "${header}" NAME)

        # skip these headers
        if (header_filename MATCHES "all.hpp|concept.hpp")
            continue ()
        endif ()

        message (STATUS "header '${header}' is missing a test")
    endforeach ()
endfunction ()

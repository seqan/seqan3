# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

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

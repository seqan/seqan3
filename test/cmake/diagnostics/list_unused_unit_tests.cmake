# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

set (seqan3_test_targets
     ""
     CACHE STRING "" FORCE)

# Add the `target` to the list of used test targets. This effectively marks the `target` as a used test.
function (collect_used_test target)
    set (seqan3_test_targets
         "${seqan3_test_targets};${target}"
         CACHE STRING "" FORCE)
endfunction ()

# Glob all test files (e.g. *.cpp files) and compare them to the list of used test targets.
function (list_unused_unit_tests)
    file (GLOB_RECURSE test_source_glob_list *.cpp)
    set (test_source_declared_list "")

    # get the source location of each "used" test target and collect it.
    foreach (test_target ${seqan3_test_targets})
        get_target_property (sources "${test_target}" SOURCES)
        get_target_property (source_dir "${test_target}" SOURCE_DIR)

        list (APPEND test_source_declared_list "${source_dir}/${sources}")
    endforeach ()

    # create the difference between test_source_glob_list set and test_source_declared_list set.
    list (REMOVE_ITEM test_source_glob_list ${test_source_declared_list})

    # list all unused tests
    foreach (test_source ${test_source_glob_list})
        message (WARNING "'${test_source}' test exists, but was not defined by seqan3_test(...)")
    endforeach ()
endfunction ()

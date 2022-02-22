# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
# -----------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.10)

set (seqan3_test_snippets
     ""
     CACHE STRING "" FORCE)

include ("${CMAKE_CURRENT_LIST_DIR}/../seqan3_path_longest_stem.cmake")
include ("${CMAKE_CURRENT_LIST_DIR}/../seqan3_test_files.cmake")

# Add the `target` to the list of used test targets. This effectively marks the `target` as a used test.
function (collect_used_snippet target)
    set (seqan3_test_snippets
         "${seqan3_test_snippets};${target}"
         CACHE STRING "" FORCE)
endfunction ()

# Glob all snippet output files (e.g. *.out and *.err files) and compare them to the list of used snippet outputs.
function (list_unused_snippets snippet_base_path)
    seqan3_test_files (test_snippet_output_glob_list "${snippet_base_path}" "*.out;*.err")
    set (test_snippet_output_list "")

    # get the source location of each "used" test target and collect it.
    foreach (test_target ${seqan3_test_snippets})
        # e.g. /seqan3/test/snippet/../../doc/tutorial/pairwise_alignment/configurations.cpp
        get_target_property (source "${test_target}" SOURCES)
        # e.g. configurations
        seqan3_path_longest_stem (source_wle "${source}")
        # e.g. /seqan3/test/snippet/../../doc/tutorial/pairwise_alignment
        get_filename_component (source_dir "${source}" DIRECTORY)
        # e.g. ../../doc/tutorial/pairwise_alignment
        file (RELATIVE_PATH source_relative_stem "${snippet_base_path}" "${source_dir}/${source_wle}")

        # test_snippet_output_list adds all potential cout / cerr files even if they don't really exist
        # This list will be subtracted from the real list of files, so it can contain "more" without changing the
        # result.
        list (APPEND test_snippet_output_list "${source_relative_stem}.out")
        list (APPEND test_snippet_output_list "${source_relative_stem}.err")
    endforeach ()

    # create the difference between test_snippet_output_glob_list set and test_snippet_output_list set.
    list (REMOVE_ITEM test_snippet_output_glob_list ${test_snippet_output_list})

    # list all unused tests
    foreach (test_snippet_output ${test_snippet_output_glob_list})
        message (AUTHOR_WARNING "'${snippet_base_path}/${test_snippet_output}' snippet output exists, "
                                "but the corresponding .cpp file is missing!")
    endforeach ()
endfunction ()

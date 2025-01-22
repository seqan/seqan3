# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

# A compatible function for cmake < 3.20 that basically returns `cmake_path (GET <filename> STEM LAST_ONLY <out_var>)`
function (seqan3_path_longest_stem out_var filename)
    cmake_path (GET filename STEM LAST_ONLY result)
    set ("${out_var}"
         "${result}"
         PARENT_SCOPE) # out-var
endfunction ()

# ======
# TESTS
# ======

seqan3_path_longest_stem (seqan3_cmake_test_path "/a/b/c/")
if (NOT seqan3_cmake_test_path STREQUAL "")
    message (FATAL_ERROR "internal error: '${seqan3_cmake_test_path}' vs '', "
                         "seqan3_path_longest_stem produces wrong result")
endif ()

seqan3_path_longest_stem (seqan3_cmake_test_path "/a/b/c/hello")
if (NOT seqan3_cmake_test_path STREQUAL "hello")
    message (FATAL_ERROR "internal error: '${seqan3_cmake_test_path}' vs 'hello', "
                         "seqan3_path_longest_stem produces wrong result")
endif ()

seqan3_path_longest_stem (seqan3_cmake_test_path "/a/b/c/hello.cpp")
if (NOT seqan3_cmake_test_path STREQUAL "hello")
    message (FATAL_ERROR "internal error: '${seqan3_cmake_test_path}' vs 'hello', "
                         "seqan3_path_longest_stem produces wrong result")
endif ()

seqan3_path_longest_stem (seqan3_cmake_test_path "/a/b/c/hello.world.cpp")
if (NOT seqan3_cmake_test_path STREQUAL "hello.world")
    message (FATAL_ERROR "internal error: '${seqan3_cmake_test_path}' vs 'hello.world', "
                         "seqan3_path_longest_stem produces wrong result")
endif ()

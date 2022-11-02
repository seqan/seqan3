# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
# -----------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.16)

# Exposes the google-benchmark target `benchmark` and `benchmark_main`.
# CMake 3.24: https://cmake.org/cmake/help/latest/module/FetchContent.html#variable:FETCHCONTENT_TRY_FIND_PACKAGE_MODE
macro (seqan3_require_benchmark)
    enable_testing ()

    set (benchmark_version "1.7.0")
    set (gbenchmark_git_tag "v${benchmark_version}")

    find_package (benchmark ${benchmark_version} EXACT QUIET)

    if (NOT benchmark_FOUND)
        message (STATUS "Fetching Google Benchmark ${benchmark_version}")

        include (FetchContent)
        FetchContent_Declare (
            gbenchmark_fetch_content
            GIT_REPOSITORY "https://github.com/google/benchmark.git"
            GIT_TAG "${gbenchmark_git_tag}")
        option (BENCHMARK_ENABLE_TESTING "" OFF)
        option (BENCHMARK_ENABLE_WERROR "" OFF) # Does not apply to Debug builds.
        option (BENCHMARK_ENABLE_INSTALL "" OFF)
        FetchContent_MakeAvailable (gbenchmark_fetch_content)
    else ()
        message (STATUS "Found Google Benchmark ${benchmark_version}")
    endif ()

    # NOTE: google benchmark's CMakeLists.txt already defines Shlwapi
    if (NOT TARGET gbenchmark_build)
        add_custom_target (gbenchmark_build DEPENDS benchmark_main benchmark)
    endif ()

endmacro ()

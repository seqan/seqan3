# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
# -----------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.16)

# Exposes the google-benchmark target `gbenchmark`.
macro (seqan3_require_benchmark)
    enable_testing ()

    set (gbenchmark_git_tag "v1.6.1")

    message (STATUS "Fetch Google Benchmark:")

    include (FetchContent)
    FetchContent_Declare (
        gbenchmark_fetch_content
        GIT_REPOSITORY "https://github.com/google/benchmark.git"
        GIT_TAG "${gbenchmark_git_tag}")
    option (BENCHMARK_ENABLE_TESTING "" OFF)
    option (BENCHMARK_ENABLE_WERROR "" OFF) # Does not apply to Debug builds.
    FetchContent_MakeAvailable (gbenchmark_fetch_content)

    # NOTE: google benchmark's CMakeLists.txt already defines Shlwapi
    add_library (gbenchmark ALIAS benchmark_main)

    if (NOT TARGET gbenchmark_build)
        add_custom_target (gbenchmark_build DEPENDS gbenchmark)
        target_compile_options ("benchmark" PUBLIC "-w")
        target_compile_options ("benchmark_main" PUBLIC "-w")
    endif ()

endmacro ()

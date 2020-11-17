# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
# -----------------------------------------------------------------------------------------------------

# Exposes the google-benchmark target `gbenchmark`.
macro (seqan3_require_benchmark_old gbenchmark_git_tag)
    set (gbenchmark_project_args ${SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS})
    list (APPEND gbenchmark_project_args "-DBENCHMARK_ENABLE_TESTING=false")
    # google-benchmarks suggest to use LTO (link-time optimisation), but we don't really need that because we are a
    # header only library. This option might be still interesting for external libraries benchmarks.
    # see https://github.com/google/benchmark#debug-vs-release
    # list (APPEND gbenchmark_project_args "-DBENCHMARK_ENABLE_LTO=true")

    # force that libraries are installed to `lib/`, because GNUInstallDirs might install it into `lib64/`
    list (APPEND gbenchmark_project_args "-DCMAKE_INSTALL_LIBDIR=${PROJECT_BINARY_DIR}/lib/")

    set (gbenchmark_path "${PROJECT_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}benchmark${CMAKE_STATIC_LIBRARY_SUFFIX}")

    include (ExternalProject)
    ExternalProject_Add (
        gbenchmark_project
        PREFIX gbenchmark_project
        GIT_REPOSITORY "https://github.com/google/benchmark.git"
        GIT_TAG "${gbenchmark_git_tag}"
        SOURCE_DIR "${SEQAN3_BENCHMARK_CLONE_DIR}"
        CMAKE_ARGS "${gbenchmark_project_args}"
        BUILD_BYPRODUCTS "${gbenchmark_path}"
        UPDATE_DISCONNECTED ${SEQAN3_TEST_BUILD_OFFLINE}
    )
    unset (gbenchmark_project_args)

    add_library (gbenchmark STATIC IMPORTED)
    add_dependencies (gbenchmark gbenchmark_project)
    set_target_properties (gbenchmark PROPERTIES IMPORTED_LOCATION "${gbenchmark_path}")
    set_property (TARGET gbenchmark APPEND PROPERTY INTERFACE_LINK_LIBRARIES "pthread")

    # NOTE: google benchmarks needs Shlwapi (Shell Lightweight Utility Functions) on windows
    # see https://msdn.microsoft.com/en-us/library/windows/desktop/bb759844(v=vs.85).aspx
    # see https://github.com/google/benchmark/blob/c614dfc0d4eadcd19b188ff9c7e226c138f894a1/README.md#platform-specific-libraries
    if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
        target_link_libraries (gbenchmark INTERFACE "Shlwapi")
    endif()

    unset (gbenchmark_path)
endmacro ()

macro (seqan3_require_benchmark)
    enable_testing ()

    set (gbenchmark_git_tag "v1.5.0")

    if (NOT CMAKE_VERSION VERSION_LESS 3.14)
        # NOTE: setting the standard prevents an ICE involving passing -fconcepts through CMAKE_CXX_FLAGS
        set (CMAKE_CXX_STANDARD 17)
        message (STATUS "Fetch google benchmark:")

        include (FetchContent)
        FetchContent_Declare (
            gbenchmark_fetch_content
            GIT_REPOSITORY "https://github.com/google/benchmark.git"
            GIT_TAG "${gbenchmark_git_tag}"
        )
        option (BENCHMARK_ENABLE_TESTING "" OFF)
        FetchContent_MakeAvailable(gbenchmark_fetch_content)

        # NOTE: google benchmark's CMakeLists.txt already defines Shlwapi
        add_library (gbenchmark ALIAS benchmark_main)
        unset (CMAKE_CXX_STANDARD)
    else ()
        message (STATUS "Use google benchmark as external project:")

        seqan3_require_benchmark_old ("${gtest_git_tag}")
    endif ()

    add_custom_target (gbenchmark_build DEPENDS gbenchmark)
endmacro ()

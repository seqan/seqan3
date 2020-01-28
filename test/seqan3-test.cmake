# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
# -----------------------------------------------------------------------------------------------------

# This file provides functionality common to the different test modules used by
# SeqAn3. To build tests, run cmake on one of the sub-folders in this directory
# which contain a CMakeLists.txt.

cmake_minimum_required (VERSION 3.7)

# require SeqAn3 package
find_package (SeqAn3 REQUIRED
              HINTS ${CMAKE_CURRENT_LIST_DIR}/../build_system)

include (CheckCXXSourceCompiles)
include (FindPackageHandleStandardArgs)
include (FindPackageMessage)

option (SEQAN3_TEST_BUILD_OFFLINE "Skip the update step of external projects." OFF)

# ----------------------------------------------------------------------------
# Paths to folders.
# ----------------------------------------------------------------------------

find_path (SEQAN3_TEST_INCLUDE_DIR NAMES seqan3/test/tmp_filename.hpp HINTS "${CMAKE_CURRENT_LIST_DIR}/include/")

set (SEQAN3_BENCHMARK_CLONE_DIR "${PROJECT_BINARY_DIR}/vendor/benchmark")
set (SEQAN3_TEST_CLONE_DIR "${PROJECT_BINARY_DIR}/vendor/googletest")

# needed for add_library (seqan3::test::* INTERFACE IMPORTED)
# see cmake bug https://gitlab.kitware.com/cmake/cmake/issues/15052
file(MAKE_DIRECTORY ${SEQAN3_BENCHMARK_CLONE_DIR}/include/)
file(MAKE_DIRECTORY ${SEQAN3_TEST_CLONE_DIR}/googletest/include/)

# ----------------------------------------------------------------------------
# Interface targets for the different test modules in seqan3.
# ----------------------------------------------------------------------------

# seqan3::test exposes a base set of required flags, includes, definitions and
# libraries which are in common for **all** seqan3 tests
add_library (seqan3_test INTERFACE)
target_compile_options (seqan3_test INTERFACE "-pedantic"  "-Wall" "-Wextra" "-Werror")
target_link_libraries (seqan3_test INTERFACE "seqan3::seqan3" "pthread")
target_include_directories (seqan3_test INTERFACE "${SEQAN3_TEST_INCLUDE_DIR}")
add_library (seqan3::test ALIAS seqan3_test)

# seqan3::test::performance specifies required flags, includes and libraries
# needed for performance test cases in seqan3/test/performance
add_library (seqan3_test_performance INTERFACE)
target_link_libraries (seqan3_test_performance INTERFACE "seqan3::test" "gbenchmark")

# NOTE: google benchmarks needs Shlwapi (Shell Lightweight Utility Functions) on windows
# see https://msdn.microsoft.com/en-us/library/windows/desktop/bb759844(v=vs.85).aspx
# see https://github.com/google/benchmark/blob/c614dfc0d4eadcd19b188ff9c7e226c138f894a1/README.md#platform-specific-libraries
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    target_link_libraries (seqan3_test_performance INTERFACE "Shlwapi")
endif()

target_include_directories (seqan3_test_performance INTERFACE "${SEQAN3_BENCHMARK_CLONE_DIR}/include/")
add_library (seqan3::test::performance ALIAS seqan3_test_performance)

# seqan3::test::unit specifies required flags, includes and libraries
# needed for unit test cases in seqan3/test/unit
add_library (seqan3_test_unit INTERFACE)
target_link_libraries (seqan3_test_unit INTERFACE "seqan3::test" "gtest_main" "gtest")
target_include_directories (seqan3_test_unit INTERFACE "${SEQAN3_TEST_CLONE_DIR}/googletest/include/")
add_library (seqan3::test::unit ALIAS seqan3_test_unit)

# seqan3::test::coverage specifies required flags, includes and libraries
# needed for coverage test cases in seqan3/test/coverage
add_library (seqan3_test_coverage INTERFACE)
target_compile_options (seqan3_test_coverage INTERFACE "--coverage" "-fprofile-arcs" "-ftest-coverage")
target_link_libraries (seqan3_test_coverage INTERFACE "seqan3::test::unit" "gcov")
add_library (seqan3::test::coverage ALIAS seqan3_test_coverage)

# seqan3::test::header specifies required flags, includes and libraries
# needed for header test cases in seqan3/test/header
add_library (seqan3_test_header INTERFACE)
target_link_libraries (seqan3_test_header INTERFACE "seqan3::test::unit")
target_link_libraries (seqan3_test_header INTERFACE "seqan3::test::performance")
add_library (seqan3::test::header ALIAS seqan3_test_header)

# ----------------------------------------------------------------------------
# Commonly shared options for external projects.
# ----------------------------------------------------------------------------

set (SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_VERBOSE_MAKEFILE=${CMAKE_VERBOSE_MAKEFILE}")

# ----------------------------------------------------------------------------
# Commonly used macros for the different test modules in seqan3.
# ----------------------------------------------------------------------------

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

macro (seqan3_test_files VAR test_base_path_ extension_wildcards)
    # test_base_path is /home/.../seqan3/test/
    get_filename_component(test_base_path "${test_base_path_}" ABSOLUTE)
    file (RELATIVE_PATH test_base_path_relative "${CMAKE_SOURCE_DIR}" "${test_base_path}")
    # ./ is a hack to deal with empty test_base_path_relative
    set (test_base_path_relative "./${test_base_path_relative}")
    # collect all cpp files
    set (${VAR} "")
    foreach (extension_wildcard ${extension_wildcards})
        file (GLOB_RECURSE test_files RELATIVE "${test_base_path}"
        "${test_base_path_relative}/${extension_wildcard}")
        list (APPEND ${VAR} ${test_files})
    endforeach ()

    unset (test_base_path)
    unset (test_base_path_relative)
endmacro ()

macro (seqan3_require_ccache)
    find_program (CCACHE_PROGRAM ccache)
    find_package_message (CCACHE_PROGRAM_PRE "Finding program ccache" "[${CCACHE_PROGRAM}]")

    if (NOT CCACHE_PROGRAM)
        find_package_message (CCACHE_PROGRAM "Finding program ccache - Failed" "[${CCACHE_PROGRAM}]")
    else ()
        find_package_message (CCACHE_PROGRAM "Finding program ccache - Success" "[${CCACHE_PROGRAM}]")
        # New option since cmake >= 3.4:
        # https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_COMPILER_LAUNCHER.html
        set (CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")

        # use ccache in external cmake projects
        list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_CXX_COMPILER_LAUNCHER=${CMAKE_CXX_COMPILER_LAUNCHER}")
    endif ()
    unset (CCACHE_PROGRAM)
endmacro ()

macro (seqan3_require_benchmark)
    enable_testing ()

    set (gbenchmark_project_args ${SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS})
    list (APPEND gbenchmark_project_args "-DBENCHMARK_ENABLE_TESTING=false")
    # list (APPEND gbenchmark_project_args "-DBENCHMARK_ENABLE_LTO=true")

    # force that libraries are installed to `lib/`, because GNUInstallDirs might install it into `lib64/`
    list (APPEND gbenchmark_project_args "-DCMAKE_INSTALL_LIBDIR=${PROJECT_BINARY_DIR}/lib/")

    set (gbenchmark_path "${PROJECT_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}benchmark${CMAKE_STATIC_LIBRARY_SUFFIX}")

    include (ExternalProject)
    ExternalProject_Add (
        gbenchmark_project
        PREFIX gbenchmark_project
        GIT_REPOSITORY "https://github.com/google/benchmark.git"
        GIT_TAG "v1.5.0"
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

    unset (gbenchmark_path)
endmacro ()

macro (seqan3_require_test)
    enable_testing ()

    set (gtest_project_args ${SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS})
    list (APPEND gtest_project_args "-DBUILD_GMOCK=0")

    # force that libraries are installed to `lib/`, because GNUInstallDirs might install it into `lib64/`
    list (APPEND gtest_project_args "-DCMAKE_INSTALL_LIBDIR=${PROJECT_BINARY_DIR}/lib/")

    # google sets CMAKE_DEBUG_POSTFIX = "d"
    set (gtest_main_path "${PROJECT_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}gtest_main${CMAKE_STATIC_LIBRARY_SUFFIX}")
    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
        set (gtest_main_path "${PROJECT_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}gtest_maind${CMAKE_STATIC_LIBRARY_SUFFIX}")
    endif ()

    set (gtest_path "${PROJECT_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}gtest${CMAKE_STATIC_LIBRARY_SUFFIX}")
    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
        set (gtest_path "${PROJECT_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}gtestd${CMAKE_STATIC_LIBRARY_SUFFIX}")
    endif ()

    include (ExternalProject)
    ExternalProject_Add (
        gtest_project
        PREFIX gtest_project
        GIT_REPOSITORY "https://github.com/google/googletest.git"
        # we currently have warnings that were introduced in
        # 03867b5389516a0f185af52672cf5472fa0c159c, which are still available
        # in "release-1.8.1", see https://github.com/google/googletest/issues/1419
        GIT_TAG "release-1.10.0"
        SOURCE_DIR "${SEQAN3_TEST_CLONE_DIR}"
        CMAKE_ARGS "${gtest_project_args}"
        BUILD_BYPRODUCTS "${gtest_main_path}" "${gtest_path}"
        UPDATE_DISCONNECTED ${SEQAN3_TEST_BUILD_OFFLINE}
    )
    unset (gtest_project_args)

    add_library (gtest_main STATIC IMPORTED)
    add_dependencies (gtest_main gtest_project)
    set_target_properties (gtest_main PROPERTIES IMPORTED_LOCATION "${gtest_main_path}")

    add_library (gtest STATIC IMPORTED)
    add_dependencies (gtest gtest_main)
    set_target_properties (gtest PROPERTIES IMPORTED_LOCATION "${gtest_path}")
    set_property (TARGET gtest APPEND PROPERTY INTERFACE_LINK_LIBRARIES "pthread")

    unset(gtest_main_path)
    unset(gtest_path)
endmacro ()

macro (add_subdirectories_of directory)
    file (GLOB ENTRIES
          RELATIVE ${directory}
          ${directory}/[!.]*)

    foreach (ENTRY ${ENTRIES})
        if (IS_DIRECTORY ${directory}/${ENTRY})
            if (EXISTS ${directory}/${ENTRY}/CMakeLists.txt)
                add_subdirectory (${directory}/${ENTRY} ${CMAKE_CURRENT_BINARY_DIR}/${ENTRY})
            endif ()
        endif ()
    endforeach ()
    unset (ENTRIES)
endmacro ()

macro (add_subdirectories)
    add_subdirectories_of(${CMAKE_CURRENT_SOURCE_DIR})
endmacro ()

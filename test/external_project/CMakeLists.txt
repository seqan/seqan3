# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
# -----------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.4)
project (seqan3_test_external_project CXX)

include (../seqan3-test.cmake) # for SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS, SEQAN3_VERSION
include (ExternalProject)

set (SEQAN3_ROOT "${CMAKE_CURRENT_LIST_DIR}/../../")

if (NOT ${CMAKE_VERSION} VERSION_LESS 3.14) # cmake 3.14 version is needed to install seqan3.
    include (install-seqan3.cmake)
endif ()

option (SEQAN3_EXTERNAL_PROJECT_FIND_DEBUG_MODE
        "Enable this option if you want to get a detailed list which paths were considered for find_package(...)"
        false)

if (NOT ${CMAKE_VERSION} VERSION_LESS 3.14) # cmake 3.14 version is needed to install cmake 3.4.
    set (SEQAN3_EXTERNAL_PROJECT_CMAKE_COMMAND "${CMAKE_COMMAND}")
    include (use-cmake3.4.cmake) # ensure that seqan3 can be used with cmake 3.4
else ()
    set (SEQAN3_EXTERNAL_PROJECT_CMAKE_COMMAND "${CMAKE_COMMAND}")
endif ()

# Here is an explanation of what happens:
#
# If your current CMake version is < 3.14, we
#  * test test/external_project/seqan3_submodule_add_subdirectory/CMakeLists.txt with your current CMake version
#  * test test/external_project/seqan3_submodule_find_package/CMakeLists.txt with your current CMake version
#  * DO NOT test test/external_project/seqan3_installed/CMakeLists.txt, because we need 3.14 to be able to
#    package seqan3 to execute the test.
#  * DO NOT test test/external_project/seqan3_fetch_content_zip/CMakeLists.txt, because we need 3.14 in order to use
#    fetch_content
#
# If your current CMake version is >= 3.14, we
#  * download CMake  3.4
#  * test test/external_project/seqan3_submodule_add_subdirectory/CMakeLists.txt with CMake 3.4
#  * test test/external_project/seqan3_submodule_find_package/CMakeLists.txt with CMake 3.4
#  * test test/external_project/seqan3_installed/CMakeLists.txt with CMake 3.4
#  * test test/external_project/seqan3_fetch_content_zip/CMakeLists.txt with your current CMake version

# 1) This tests test/external_project/seqan3_submodule_add_subdirectory/CMakeLists.txt
#    That means we use add_subdirectory directly on seqan3's top level CMakeLists.txt.
#    This will automatically call find_package and expose our seqan3::seqan3 target.
#    This is expected to work with CMake >= 3.4.
# (ExternalProject_Add simulates a fresh and separate invocation of cmake ../)
ExternalProject_Add (
    seqan3_submodule_add_subdirectory
    PREFIX seqan3_submodule_add_subdirectory
    SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/seqan3_submodule_add_subdirectory"
    CMAKE_COMMAND "${SEQAN3_EXTERNAL_PROJECT_CMAKE_COMMAND}"
    CMAKE_ARGS
        ${SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS}
        "-DSEQAN3_ROOT=${SEQAN3_ROOT}"
        "-DCMAKE_FIND_DEBUG_MODE=${SEQAN3_EXTERNAL_PROJECT_FIND_DEBUG_MODE}"
)

# 2) This tests test/external_project/seqan3_submodule_find_package/CMakeLists.txt
#    We have a seqan3 checkout somewhere and we point CMAKE_PREFIX_PATH to <checkout>/seqan3/build_system
#    and then use `find_package` to find `seqan3-config.cmake` which exposes our `seqan3::seqan3` target.
#    This is expected to work with CMake >= 3.4.
# (ExternalProject_Add simulates a fresh and separate invocation of cmake ../)
ExternalProject_Add (
    seqan3_submodule_find_package
    PREFIX seqan3_submodule_find_package
    SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/seqan3_submodule_find_package"
    CMAKE_COMMAND "${SEQAN3_EXTERNAL_PROJECT_CMAKE_COMMAND}"
    CMAKE_ARGS
        ${SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS}
        "-DCMAKE_PREFIX_PATH=${SEQAN3_ROOT}/build_system"
        "-DCMAKE_FIND_DEBUG_MODE=${SEQAN3_EXTERNAL_PROJECT_FIND_DEBUG_MODE}"
)

# 3) This tests test/external_project/seqan3_installed/CMakeLists.txt
#    This test assumes that seqan3 was installed by make install (e.g. system-wide).
#    This is the way most upstream packages, like debian, provide our library.
#    This test assumes that `seqan3-config.cmake` can be found by cmake in some global paths like /usr/share/cmake/.
#
#    We simulate this by using our `make package` release, e.g. the one we release under
#    https://github.com/seqan/seqan3/releases, and unzipping it to some folder and making
#    that path globally accessible by CMAKE_SYSTEM_PREFIX_PATH.
#    We need CMake >= 3.14 to be able to package seqan3, but we actually expect that this
#    test works with CMake >= 3.4.
# (ExternalProject_Add simulates a fresh and separate invocation of cmake ../)
if (NOT ${CMAKE_VERSION} VERSION_LESS 3.14) # cmake 3.14 version is needed to install seqan3.
    ExternalProject_Add (
        seqan3_installed
        PREFIX seqan3_installed
        SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/seqan3_installed"
        CMAKE_COMMAND "${SEQAN3_EXTERNAL_PROJECT_CMAKE_COMMAND}"
        CMAKE_ARGS
            ${SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS}
            "-DCMAKE_SYSTEM_PREFIX_PATH=${SEQAN3_SYSTEM_PREFIX}"
            "-DCMAKE_FIND_DEBUG_MODE=${SEQAN3_EXTERNAL_PROJECT_FIND_DEBUG_MODE}"
    )
    add_dependencies(seqan3_installed seqan3_test_prerequisite)
endif ()

# 4) This tests test/external_project/seqan3_fetch_content_zip/CMakeLists.txt
#    It uses fetch_content (a CMake 3.14 feature) to download our zip-release (e.g. zip, tar.xz) from
#    https://github.com/seqan/seqan3/releases. fetch_content will automatically download, verify, extract it.
#    The user only needs to define CMAKE_PREFIX_PATH to be able to find our `seqan3-config.cmake`.
#    Note that FetchContent is a CMake >= 3.14 feature.
#    This is expected to work with CMake >= 3.14.
# (ExternalProject_Add simulates a fresh and separate invocation of cmake ../)
if (NOT ${CMAKE_VERSION} VERSION_LESS 3.14) # cmake 3.14 version is needed to use fetch_content.
    ExternalProject_Add (
        seqan3_fetch_content_zip
        PREFIX seqan3_fetch_content_zip
        SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/seqan3_fetch_content_zip"
        CMAKE_ARGS
            ${SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS}
            "-DSEQAN3_PACKAGE_ZIP_URL=${SEQAN3_PACKAGE_ZIP_URL}"
            "-DCMAKE_FIND_DEBUG_MODE=${SEQAN3_EXTERNAL_PROJECT_FIND_DEBUG_MODE}"
    )
    add_dependencies(seqan3_fetch_content_zip seqan3_test_prerequisite)
endif ()

# 5) This test is the same as 2) but emulates the settings within the setup tutorial.
#    This test is used as snippet in the setup tutorial.
ExternalProject_Add (
    seqan3_setup_tutorial
    PREFIX seqan3_setup_tutorial
    SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/seqan3_setup_tutorial"
    CMAKE_COMMAND "${SEQAN3_EXTERNAL_PROJECT_CMAKE_COMMAND}"
    CMAKE_ARGS
        ${SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS}
        "-DCMAKE_FIND_DEBUG_MODE=${SEQAN3_EXTERNAL_PROJECT_FIND_DEBUG_MODE}"
)

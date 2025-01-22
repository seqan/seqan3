# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.20...3.31)
project (seqan3_test_external_project CXX)

include (../seqan3-test.cmake) # for SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS, SEQAN3_VERSION
include (ExternalProject)

set (SEQAN3_ROOT "${CMAKE_CURRENT_LIST_DIR}/../../")

include (install-seqan3.cmake)
include (install-sharg.cmake)

option (SEQAN3_EXTERNAL_PROJECT_FIND_DEBUG_MODE
        "Enable this option if you want to get a detailed list which paths were considered for find_package(...)" false)

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
#  * download CMake  3.5
#  * test test/external_project/seqan3_submodule_add_subdirectory/CMakeLists.txt with CMake 3.5
#  * test test/external_project/seqan3_submodule_find_package/CMakeLists.txt with CMake 3.5
#  * test test/external_project/seqan3_installed/CMakeLists.txt with CMake 3.5
#  * test test/external_project/seqan3_fetch_content_zip/CMakeLists.txt with your current CMake version

# 1) This tests test/external_project/seqan3_submodule_add_subdirectory/CMakeLists.txt
#    That means we use add_subdirectory directly on seqan3's top level CMakeLists.txt.
#    This will automatically call find_package and expose our seqan3::seqan3 target.
#    This is expected to work with CMake >= 3.5.
# (ExternalProject_Add simulates a fresh and separate invocation of cmake ../)
ExternalProject_Add (
    seqan3_submodule_add_subdirectory
    PREFIX seqan3_submodule_add_subdirectory
    SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/seqan3_submodule_add_subdirectory"
    CMAKE_ARGS ${SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS}
               "-DCMAKE_FIND_DEBUG_MODE=${SEQAN3_EXTERNAL_PROJECT_FIND_DEBUG_MODE}" #
               "-DSEQAN3_ROOT=${SEQAN3_ROOT}")

# 2) This tests test/external_project/seqan3_submodule_find_package/CMakeLists.txt
#    We have a seqan3 checkout somewhere and we point CMAKE_PREFIX_PATH to <checkout>/seqan3/cmake
#    and then use `find_package` to find `seqan3-config.cmake` which exposes our `seqan3::seqan3` target.
#    This is expected to work with CMake >= 3.5.
# (ExternalProject_Add simulates a fresh and separate invocation of cmake ../)
ExternalProject_Add (
    seqan3_submodule_find_package
    PREFIX seqan3_submodule_find_package
    SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/seqan3_submodule_find_package"
    CMAKE_ARGS ${SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS}
               "-DCMAKE_FIND_DEBUG_MODE=${SEQAN3_EXTERNAL_PROJECT_FIND_DEBUG_MODE}"
               "-DCMAKE_PREFIX_PATH=${SEQAN3_ROOT}/cmake")

# 3) This tests test/external_project/seqan3_installed/CMakeLists.txt
#    This test assumes that seqan3 was installed by make install (e.g. system-wide).
#    This is the way most upstream packages, like debian, provide our library.
#    This test assumes that `seqan3-config.cmake` can be found by cmake in some global paths like /usr/share/cmake/.
#
#    We simulate this by using our `make package` release, e.g. the one we release under
#    https://github.com/seqan/seqan3/releases, and unzipping it to some folder and making
#    that path globally accessible by CMAKE_SYSTEM_PREFIX_PATH.
#    We need CMake >= 3.14 to be able to package seqan3, but we actually expect that this
#    test works with CMake >= 3.5.
# (ExternalProject_Add simulates a fresh and separate invocation of cmake ../)
ExternalProject_Add (
    seqan3_installed
    PREFIX seqan3_installed
    SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/seqan3_installed"
    CMAKE_ARGS ${SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS}
               "-DCMAKE_FIND_DEBUG_MODE=${SEQAN3_EXTERNAL_PROJECT_FIND_DEBUG_MODE}"
               "-DCMAKE_SYSTEM_PREFIX_PATH=${SEQAN3_SYSTEM_PREFIX}")
add_dependencies (seqan3_installed seqan3_test_prerequisite)

# 4) This tests test/external_project/seqan3_fetch_content_zip/CMakeLists.txt
#    It uses fetch_content (a CMake 3.14 feature) to download our zip-release (e.g. zip, tar.xz) from
#    https://github.com/seqan/seqan3/releases. fetch_content will automatically download, verify, extract it.
#    The user only needs to define CMAKE_PREFIX_PATH to be able to find our `seqan3-config.cmake`.
#    Note that FetchContent is a CMake >= 3.14 feature.
#    This is expected to work with CMake >= 3.14.
# (ExternalProject_Add simulates a fresh and separate invocation of cmake ../)
ExternalProject_Add (
    seqan3_fetch_content_zip
    PREFIX seqan3_fetch_content_zip
    SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/seqan3_fetch_content_zip"
    CMAKE_ARGS ${SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS}
               "-DCMAKE_FIND_DEBUG_MODE=${SEQAN3_EXTERNAL_PROJECT_FIND_DEBUG_MODE}"
               "-DSEQAN3_PACKAGE_ZIP_URL=${SEQAN3_PACKAGE_ZIP_URL}")
add_dependencies (seqan3_fetch_content_zip seqan3_test_prerequisite)

# 5) This test is the same as 2) but emulates the settings within the setup tutorial.
#    This test is used as snippet in the setup tutorial.
ExternalProject_Add (
    seqan3_setup_tutorial
    PREFIX seqan3_setup_tutorial
    SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/seqan3_setup_tutorial"
    CMAKE_ARGS ${SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS}
               "-DCMAKE_FIND_DEBUG_MODE=${SEQAN3_EXTERNAL_PROJECT_FIND_DEBUG_MODE}")

# 6) This test is the same as 5), but with Sharg via checkout.
ExternalProject_Add (
    seqan3_setup_tutorial_with_sharg
    PREFIX seqan3_setup_tutorial_with_sharg
    SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/seqan3_setup_tutorial_with_sharg"
    CMAKE_ARGS ${SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS}
               "-DCMAKE_FIND_DEBUG_MODE=${SEQAN3_EXTERNAL_PROJECT_FIND_DEBUG_MODE}" #
               "-DSHARG_ROOT=${SHARG_ROOT}")
add_dependencies (seqan3_setup_tutorial_with_sharg sharg_test_prerequisite)

# 7) This test is the same as 5), but with Sharg and SeqAn3 installed.
ExternalProject_Add (
    seqan3_setup_tutorial_with_sharg_installed
    PREFIX seqan3_setup_tutorial_with_sharg_installed
    SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/seqan3_setup_tutorial_with_sharg_installed"
    CMAKE_ARGS ${SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS}
               "-DCMAKE_FIND_DEBUG_MODE=${SEQAN3_EXTERNAL_PROJECT_FIND_DEBUG_MODE}"
               "-DCMAKE_SYSTEM_PREFIX_PATH=${SEQAN3_SYSTEM_PREFIX}")
add_dependencies (seqan3_setup_tutorial_with_sharg_installed seqan3_test_prerequisite sharg_test_prerequisite)

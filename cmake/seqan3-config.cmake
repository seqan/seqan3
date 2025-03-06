# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause
#
# This CMake module will try to find SeqAn and its dependencies.  You can use
# it the same way you would use any other CMake module.
#
#   find_package (SeqAn3 [REQUIRED] ...)
#
# Since this makes a difference for CMAKE, pay attention to the case
# ("SeqAn3", "SEQAN3" and "seqan3" are all valid, but other names not).
#
# SeqAn has the following platform requirements:
#
#   C++20
#   pthread
#
# SeqAn has the following optional dependencies:
#
#   ZLIB      -- zlib compression library
#   BZip2     -- libbz2 compression library
#   Cereal    -- Serialisation library
#
# If you don't wish for these to be detected (and used), you may define SEQAN3_NO_ZLIB,
# SEQAN3_NO_BZIP2, and SEQAN3_NO_CEREAL respectively.
#
# If you wish to require the presence of ZLIB or BZip2, just check for the module before
# finding SeqAn3, e.g. "find_package (ZLIB REQUIRED)" and "find_package (BZip2 REQUIRED)".
# If you wish to require the presence of CEREAL, you may define SEQAN3_CEREAL.
#
# Once the search has been performed, the following variables will be set.
#
#   SEQAN3_FOUND            -- Indicate whether SeqAn was found and requirements met.
#
#   SEQAN3_VERSION          -- The version as string, e.g. "3.0.0"
#   SEQAN3_VERSION_MAJOR    -- e.g. 3
#   SEQAN3_VERSION_MINOR    -- e.g. 0
#   SEQAN3_VERSION_PATCH    -- e.g. 0
#
#   SEQAN3_INCLUDE_DIRS     -- to be passed to include_directories ()
#   SEQAN3_LIBRARIES        -- to be passed to target_link_libraries ()
#   SEQAN3_DEFINITIONS      -- to be passed to add_definitions ()
#
# Additionally, the following [IMPORTED][IMPORTED] targets are defined:
#
#   seqan3::seqan3          -- interface target where
#                                  target_link_libraries(target seqan3::seqan3)
#                              automatically sets
#                                  target_include_directories(target $SEQAN3_INCLUDE_DIRS),
#                                  target_link_libraries(target $SEQAN3_LIBRARIES) and
#                                  target_compile_definitions(target $SEQAN3_DEFINITIONS)
#                              for a target.
#
#   [IMPORTED]: https://cmake.org/cmake/help/v3.10/prop_tgt/IMPORTED.html#prop_tgt:IMPORTED
#
# ============================================================================

# ----------------------------------------------------------------------------
# Set initial variables
# ----------------------------------------------------------------------------

# make output globally quiet if required by find_package, this effects cmake functions like `check_*`
set (CMAKE_REQUIRED_QUIET_SAVE ${CMAKE_REQUIRED_QUIET})
set (CMAKE_REQUIRED_QUIET 1)

# ----------------------------------------------------------------------------
# Greeter
# ----------------------------------------------------------------------------

string (ASCII 27 Esc)
set (ColourBold "${Esc}[1m")
set (ColourReset "${Esc}[m")

if (NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
    message (STATUS "${ColourBold}Finding SeqAn3 and checking requirements:${ColourReset}")
endif ()

# ----------------------------------------------------------------------------
# Includes
# ----------------------------------------------------------------------------

include (CheckIncludeFileCXX)
include (CheckCXXSourceCompiles)
include (FindPackageHandleStandardArgs)

# ----------------------------------------------------------------------------
# Pretty printing and error handling
# ----------------------------------------------------------------------------

macro (seqan3_config_print text)
    if (NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
        message (STATUS "  ${text}")
    endif ()
endmacro ()

macro (seqan3_config_error text)
    if (${CMAKE_FIND_PACKAGE_NAME}_FIND_REQUIRED)
        message (FATAL_ERROR ${text})
    else ()
        if (NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
            message (WARNING ${text})
        endif ()
        return ()
    endif ()
endmacro ()

# ----------------------------------------------------------------------------
# CPM
# ----------------------------------------------------------------------------

# This will be true for git clones and source packages, but not for installed packages.
if (EXISTS "${CMAKE_CURRENT_LIST_DIR}/CPM.cmake")
    set (SEQAN3_HAS_CPM TRUE)
else ()
    set (SEQAN3_HAS_CPM FALSE)
endif ()

if (SEQAN3_HAS_CPM)
    set (CPM_INDENT "  CMake Package Manager CPM: ")
    include ("${CMAKE_CURRENT_LIST_DIR}/CPM.cmake")
    CPMUsePackageLock ("${CMAKE_CURRENT_LIST_DIR}/package-lock.cmake")
endif ()

# ----------------------------------------------------------------------------
# Find SeqAn3 include path
# ----------------------------------------------------------------------------

# Note that seqan3-config.cmake can be standalone and thus SEQAN3_CLONE_DIR might be empty.
# * `SEQAN3_INCLUDE_DIR` was already found in seqan3-config-version.cmake
if (SEQAN3_INCLUDE_DIR)
    seqan3_config_print ("SeqAn3 include dir found:   ${SEQAN3_INCLUDE_DIR}")
else ()
    seqan3_config_error ("SeqAn3 include directory could not be found (SEQAN3_INCLUDE_DIR: '${SEQAN3_INCLUDE_DIR}')")
endif ()

# ----------------------------------------------------------------------------
# Force-(de)activate optional dependencies
# ----------------------------------------------------------------------------

# https://cmake.org/cmake/help/latest/variable/CMAKE_DISABLE_FIND_PACKAGE_PackageName.html
# https://cmake.org/cmake/help/latest/variable/CMAKE_REQUIRE_FIND_PACKAGE_PackageName.html

## Example for deactivating
# cmake <path> -DCMAKE_DISABLE_FIND_PACKAGE_ZLIB=TRUE
#              -DCMAKE_DISABLE_FIND_PACKAGE_BZip2=TRUE
#              -DCMAKE_DISABLE_FIND_PACKAGE_cereal=TRUE

## Example for requiring
# cmake <path> -DCMAKE_REQUIRE_FIND_PACKAGE_ZLIB=TRUE
#              -DCMAKE_REQUIRE_FIND_PACKAGE_BZip2=TRUE
#              -DCMAKE_REQUIRE_FIND_PACKAGE_cereal=TRUE

# ----------------------------------------------------------------------------
# Check supported compilers
# ----------------------------------------------------------------------------

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 12)
    message (FATAL_ERROR "GCC < 12 is not supported. The detected compiler version is ${CMAKE_CXX_COMPILER_VERSION}.")
endif ()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 17)
    message (FATAL_ERROR "Clang < 17 is not supported. The detected compiler version is ${CMAKE_CXX_COMPILER_VERSION}.")
endif ()

# ----------------------------------------------------------------------------
# thread support (pthread, windows threads)
# ----------------------------------------------------------------------------

set (THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package (Threads QUIET)

if (TARGET Threads::Threads)
    list (APPEND SEQAN3_LIBRARIES Threads::Threads)
    if ("${CMAKE_THREAD_LIBS_INIT}" STREQUAL "")
        seqan3_config_print ("Thread support:             builtin.")
    else ()
        seqan3_config_print ("Thread support:             via ${CMAKE_THREAD_LIBS_INIT}")
    endif ()
else ()
    seqan3_config_print ("Thread support:             not found.")
endif ()

# ----------------------------------------------------------------------------
# Cereal dependency
# ----------------------------------------------------------------------------

if (SEQAN3_HAS_CPM AND NOT CMAKE_DISABLE_FIND_PACKAGE_cereal)
    CPMGetPackage (cereal)
else ()
    find_package (cereal CONFIG QUIET)
endif ()

if (TARGET cereal::cereal)
    list (APPEND SEQAN3_LIBRARIES cereal::cereal)
    seqan3_config_print ("Optional dependency:        Cereal found.")
else ()
    set (SEQAN3_DEFINITIONS ${SEQAN3_DEFINITIONS} "-DSEQAN3_HAS_CEREAL=0")
    seqan3_config_print ("Optional dependency:        Cereal not found.")
endif ()

# ----------------------------------------------------------------------------
# ZLIB dependency
# ----------------------------------------------------------------------------

find_package (ZLIB QUIET)

if (TARGET ZLIB::ZLIB)
    list (APPEND SEQAN3_LIBRARIES ZLIB::ZLIB)
    seqan3_config_print ("Optional dependency:        ZLIB-${ZLIB_VERSION_STRING} found.")
else ()
    set (SEQAN3_DEFINITIONS ${SEQAN3_DEFINITIONS} "-DSEQAN3_HAS_ZLIB=0")
    seqan3_config_print ("Optional dependency:        ZLIB not found.")
endif ()

# ----------------------------------------------------------------------------
# BZip2 dependency
# ----------------------------------------------------------------------------

find_package (BZip2 QUIET)

if (TARGET ZLIB::ZLIB AND TARGET BZip2::BZip2)
    list (APPEND SEQAN3_LIBRARIES BZip2::BZip2)
    seqan3_config_print ("Optional dependency:        BZip2-${BZIP2_VERSION_STRING} found.")
else ()
    set (SEQAN3_DEFINITIONS ${SEQAN3_DEFINITIONS} "-DSEQAN3_HAS_BZIP2=0")
    seqan3_config_print ("Optional dependency:        BZip2 not found.")
endif ()

if (NOT TARGET ZLIB::ZLIB AND TARGET BZip2::BZip2)
    message (AUTHOR_WARNING "BZip2 was found but ZLIB was not found. BZip2 requires ZLIB.")
endif ()

# ----------------------------------------------------------------------------
# System dependencies
# ----------------------------------------------------------------------------

# librt
find_library (SEQAN3_RT_LIB rt)
if (SEQAN3_RT_LIB)
    list (APPEND SEQAN3_LIBRARIES ${SEQAN3_RT_LIB})
endif ()

# libexecinfo -- implicit
find_package (Backtrace QUIET)
if (TARGET Backtrace::Backtrace)
    list (APPEND SEQAN3_LIBRARIES Backtrace::Backtrace)
    seqan3_config_print ("Optional dependency:        libexecinfo found.")
else ()
    seqan3_config_print ("Optional dependency:        libexecinfo not found.")
endif ()

# ----------------------------------------------------------------------------
# Perform compilability test of platform.hpp (tests some requirements)
# ----------------------------------------------------------------------------

# cmake-format: off
# Note: With CMake >= 3.25, the file WRITE can be removed, the second and third line in try_compile can be replaced by
# SOURCE_FROM_CONTENT "platform_test.cpp" "#include <seqan3/core/platform.hpp>\nint main() {}"
file (WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/platform_test.cpp"
            "#include <seqan3/core/platform.hpp>\nint main() {}")

try_compile (SEQAN3_PLATFORM_TEST
             ${CMAKE_BINARY_DIR}
             ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/platform_test.cpp
             CMAKE_FLAGS "-DINCLUDE_DIRECTORIES:STRING=${SEQAN3_INCLUDE_DIR}"
             COMPILE_DEFINITIONS ${SEQAN3_DEFINITIONS}
             CXX_STANDARD 23
             CXX_STANDARD_REQUIRED ON
             CXX_EXTENSIONS OFF
             OUTPUT_VARIABLE SEQAN3_PLATFORM_TEST_OUTPUT)
# cmake-format: on

if (SEQAN3_PLATFORM_TEST)
    seqan3_config_print ("SeqAn3 platform.hpp build:  passed.")
else ()
    seqan3_config_error ("SeqAn3 platform.hpp build:  failed!\n\
                        ${SEQAN3_PLATFORM_TEST_OUTPUT}")
endif ()

# ----------------------------------------------------------------------------
# Finish find_package call
# ----------------------------------------------------------------------------

if (CMAKE_FIND_PACKAGE_NAME)
    find_package_handle_standard_args (${CMAKE_FIND_PACKAGE_NAME} REQUIRED_VARS SEQAN3_INCLUDE_DIR)

    # Set SEQAN3_* variables with the content of ${CMAKE_FIND_PACKAGE_NAME}_(FOUND|...|VERSION)
    # This needs to be done, because `find_package(SeqAn3)` might be called in any case-sensitive way and we want to
    # guarantee that SEQAN3_* are always set.
    foreach (package_var
             FOUND
             DIR
             ROOT
             CONFIG
             VERSION
             VERSION_MAJOR
             VERSION_MINOR
             VERSION_PATCH
             VERSION_TWEAK
             VERSION_COUNT)
        set (SEQAN3_${package_var} "${${CMAKE_FIND_PACKAGE_NAME}_${package_var}}")
    endforeach ()
else ()
    set (SEQAN3_VERSION "${PACKAGE_VERSION}")
endif ()

# propagate SEQAN3_INCLUDE_DIR into SEQAN3_INCLUDE_DIRS
set (SEQAN3_INCLUDE_DIRS ${SEQAN3_INCLUDE_DIR})

# ----------------------------------------------------------------------------
# Export targets
# ----------------------------------------------------------------------------

if (NOT TARGET seqan3::seqan3)
    add_library (seqan3_seqan3 INTERFACE)
    target_compile_definitions (seqan3_seqan3 INTERFACE ${SEQAN3_DEFINITIONS})
    target_compile_features (seqan3_seqan3 INTERFACE cxx_std_23)
    target_link_libraries (seqan3_seqan3 INTERFACE "${SEQAN3_LIBRARIES}")
    # include seqan3/include/ as -I, because seqan3 should never produce warnings.
    target_include_directories (seqan3_seqan3 INTERFACE "${SEQAN3_INCLUDE_DIR}")
    add_library (seqan3::seqan3 ALIAS seqan3_seqan3)
endif ()

set (CMAKE_REQUIRED_QUIET ${CMAKE_REQUIRED_QUIET_SAVE})

if (SEQAN3_FIND_DEBUG)
    message ("Result for ${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt")
    message ("")
    message ("  CMAKE_BUILD_TYPE            ${CMAKE_BUILD_TYPE}")
    message ("  CMAKE_SOURCE_DIR            ${CMAKE_SOURCE_DIR}")
    message ("  CMAKE_INCLUDE_PATH          ${CMAKE_INCLUDE_PATH}")
    message ("  SEQAN3_INCLUDE_DIR          ${SEQAN3_INCLUDE_DIR}")
    message ("")
    message ("  ${CMAKE_FIND_PACKAGE_NAME}_FOUND                ${${CMAKE_FIND_PACKAGE_NAME}_FOUND}")
    message ("  SEQAN3_HAS_ZLIB             ${ZLIB_FOUND}")
    message ("  SEQAN3_HAS_BZIP2            ${BZIP2_FOUND}")
    message ("")
    message ("  SEQAN3_INCLUDE_DIRS         ${SEQAN3_INCLUDE_DIRS}")
    message ("  SEQAN3_LIBRARIES            ${SEQAN3_LIBRARIES}")
    message ("  SEQAN3_DEFINITIONS          ${SEQAN3_DEFINITIONS}")
    message ("")
    message ("  SEQAN3_VERSION              ${SEQAN3_VERSION}")
    message ("  SEQAN3_VERSION_MAJOR        ${SEQAN3_VERSION_MAJOR}")
    message ("  SEQAN3_VERSION_MINOR        ${SEQAN3_VERSION_MINOR}")
    message ("  SEQAN3_VERSION_PATCH        ${SEQAN3_VERSION_PATCH}")
endif ()

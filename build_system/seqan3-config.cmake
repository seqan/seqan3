# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
# -----------------------------------------------------------------------------------------------------
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
# SeqAn requires the following libraries:
#
#   SDSL      -- the succinct data structure library
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
#   SEQAN3_CXX_FLAGS        -- to be added to CMAKE_CXX_FLAGS
#
# Additionally, the following [IMPORTED][IMPORTED] targets are defined:
#
#   seqan3::seqan3          -- interface target where
#                                  target_link_libraries(target seqan3::seqan3)
#                              automatically sets
#                                  target_include_directories(target $SEQAN3_INCLUDE_DIRS),
#                                  target_link_libraries(target $SEQAN3_LIBRARIES),
#                                  target_compile_definitions(target $SEQAN3_DEFINITIONS) and
#                                  target_compile_options(target $SEQAN3_CXX_FLAGS)
#                              for a target.
#
#   [IMPORTED]: https://cmake.org/cmake/help/v3.10/prop_tgt/IMPORTED.html#prop_tgt:IMPORTED
#
# ============================================================================

cmake_minimum_required (VERSION 3.4...3.12)

# ----------------------------------------------------------------------------
# Set initial variables
# ----------------------------------------------------------------------------

# make output globally quiet if required by find_package, this effects cmake functions like `check_*`
set (CMAKE_REQUIRED_QUIET_SAVE ${CMAKE_REQUIRED_QUIET})
set (CMAKE_REQUIRED_QUIET ${${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY})

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
# Find SeqAn3 include path
# ----------------------------------------------------------------------------

# Note that seqan3-config.cmake can be standalone and thus SEQAN3_CLONE_DIR might be empty.
# * `SEQAN3_CLONE_DIR` was already found in seqan3-config-version.cmake
# * `SEQAN3_INCLUDE_DIR` was already found in seqan3-config-version.cmake
find_path (SEQAN3_SUBMODULES_DIR
           NAMES submodules/sdsl-lite
           HINTS "${SEQAN3_CLONE_DIR}" "${SEQAN3_INCLUDE_DIR}/seqan3")

if (SEQAN3_INCLUDE_DIR)
    seqan3_config_print ("SeqAn3 include dir found:   ${SEQAN3_INCLUDE_DIR}")
else ()
    seqan3_config_error ("SeqAn3 include directory could not be found (SEQAN3_INCLUDE_DIR: '${SEQAN3_INCLUDE_DIR}')")
endif ()

# ----------------------------------------------------------------------------
# Detect if we are a clone of repository and if yes auto-add submodules
# ----------------------------------------------------------------------------

if (SEQAN3_CLONE_DIR)
    seqan3_config_print ("Detected as running from a repository checkout…")
endif ()

if (SEQAN3_SUBMODULES_DIR)
    file (GLOB submodules ${SEQAN3_SUBMODULES_DIR}/submodules/*/include)
    foreach (submodule ${submodules})
        if (IS_DIRECTORY ${submodule})
            seqan3_config_print ("  …adding submodule include:  ${submodule}")
            set (SEQAN3_DEPENDENCY_INCLUDE_DIRS ${submodule} ${SEQAN3_DEPENDENCY_INCLUDE_DIRS})
        endif ()
    endforeach ()
endif ()

# ----------------------------------------------------------------------------
# Options for CheckCXXSourceCompiles
# ----------------------------------------------------------------------------

# deactivate messages in check_*
set (CMAKE_REQUIRED_QUIET 1)
# use global variables in Check* calls
set (CMAKE_REQUIRED_INCLUDES ${CMAKE_INCLUDE_PATH} ${SEQAN3_INCLUDE_DIR} ${SEQAN3_DEPENDENCY_INCLUDE_DIRS})
set (CMAKE_REQUIRED_FLAGS ${CMAKE_CXX_FLAGS})

# ----------------------------------------------------------------------------
# Force-deactivate optional dependencies
# ----------------------------------------------------------------------------

# Cereal is auto-detected by default, i.e. used if found, not used if not found.
# You can optionally set a hard requirement so a build fails without cereal,
# or you can force-disable cereal even if present on the system.
option (SEQAN3_CEREAL "Require cereal and fail if not present." OFF)
option (SEQAN3_NO_CEREAL "Don't use cereal, even if present." OFF)

if (SEQAN3_CEREAL AND SEQAN3_NO_CEREAL)
    # this is always a user error, therefore we always error-out, even if SeqAn is not required
    message (FATAL_ERROR "You may not specify SEQAN3_CEREAL and SEQAN3_NO_CEREAL at the same time.\n\
                          You can specify neither (use auto-detection), or specify either to force on/off.")
    return ()
endif ()

if (SEQAN3_CEREAL)
    set (SEQAN3_DEFINITIONS ${SEQAN3_DEFINITIONS} "-DSEQAN3_WITH_CEREAL=1")
elseif (SEQAN3_NO_CEREAL)
    set (SEQAN3_DEFINITIONS ${SEQAN3_DEFINITIONS} "-DSEQAN3_WITH_CEREAL=0")
endif ()

# These two are "opt-in", because detected by CMake
# If you want to force-require these, just do find_package (zlib REQUIRED) before find_package (seqan3)
option (SEQAN3_NO_ZLIB "Don't use ZLIB, even if present." OFF)
option (SEQAN3_NO_BZIP2 "Don't use BZip2, even if present." OFF)

# ----------------------------------------------------------------------------
# Check supported compilers
# ----------------------------------------------------------------------------

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 10)
    message (FATAL_ERROR "GCC < 10 is not supported. The detected compiler version is ${CMAKE_CXX_COMPILER_VERSION}.")
endif ()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 17)
    message (FATAL_ERROR "Clang < 17 is not supported. The detected compiler version is ${CMAKE_CXX_COMPILER_VERSION}.")
endif ()

# ----------------------------------------------------------------------------
# Require C++20
# ----------------------------------------------------------------------------

set (CMAKE_REQUIRED_FLAGS_SAVE ${CMAKE_REQUIRED_FLAGS})

set (CXXSTD_TEST_SOURCE
     "#if !defined (__cplusplus) || (__cplusplus < 202002L)
      #error NOCXX20
      #endif
      int main() {}")

set (SEQAN3_FEATURE_CPP20_FLAG_BUILTIN "")
set (SEQAN3_FEATURE_CPP20_FLAG_STD20 "-std=c++20")

set (SEQAN3_CPP20_FLAG "")

foreach (_FLAG BUILTIN STD20)
    set (CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS_SAVE} ${SEQAN3_FEATURE_CPP20_FLAG_${_FLAG}}")

    check_cxx_source_compiles ("${CXXSTD_TEST_SOURCE}" CPP20_FLAG_${_FLAG})

    if (CPP20_FLAG_${_FLAG})
        set (SEQAN3_CPP20_FLAG ${_FLAG})
        break ()
    endif ()
endforeach ()

set (CMAKE_REQUIRED_FLAGS ${CMAKE_REQUIRED_FLAGS_SAVE})

if (SEQAN3_CPP20_FLAG STREQUAL "BUILTIN")
    seqan3_config_print ("C++ Standard-20 support:    builtin")
elseif (SEQAN3_CPP20_FLAG)
    set (SEQAN3_CXX_FLAGS "${SEQAN3_CXX_FLAGS} ${SEQAN3_FEATURE_CPP20_FLAG_${SEQAN3_CPP20_FLAG}}")
    seqan3_config_print ("C++ Standard-20 support:    via ${SEQAN3_FEATURE_CPP20_FLAG_${SEQAN3_CPP20_FLAG}}")
else ()
    seqan3_config_error ("SeqAn3 requires C++20, but your compiler does not support it.")
endif ()

# ----------------------------------------------------------------------------
# thread support (pthread, windows threads)
# ----------------------------------------------------------------------------

set (THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package (Threads QUIET)

if (Threads_FOUND)
    set (SEQAN3_LIBRARIES ${SEQAN3_LIBRARIES} Threads::Threads)
    if ("${CMAKE_THREAD_LIBS_INIT}" STREQUAL "")
        seqan3_config_print ("Thread support:             builtin.")
    else ()
        seqan3_config_print ("Thread support:             via ${CMAKE_THREAD_LIBS_INIT}")
    endif ()
else ()
    seqan3_config_print ("Thread support:             not found.")
endif ()

# ----------------------------------------------------------------------------
# Require SDSL
# ----------------------------------------------------------------------------

check_include_file_cxx (sdsl/version.hpp _SEQAN3_HAVE_SDSL)

if (_SEQAN3_HAVE_SDSL)
    seqan3_config_print ("Required dependency:        SDSL found.")
else ()
    seqan3_config_error (
        "The SDSL library is required, but wasn't found. Get it from https://github.com/xxsds/sdsl-lite")
endif ()

# ----------------------------------------------------------------------------
# Cereal dependency is optional, but may set as required
# ----------------------------------------------------------------------------

if (NOT SEQAN3_NO_CEREAL)
    check_include_file_cxx (cereal/cereal.hpp _SEQAN3_HAVE_CEREAL)

    if (_SEQAN3_HAVE_CEREAL)
        if (SEQAN3_CEREAL)
            seqan3_config_print ("Required dependency:        Cereal found.")
        else ()
            seqan3_config_print ("Optional dependency:        Cereal found.")
        endif ()
    else ()
        if (SEQAN3_CEREAL)
            seqan3_config_error ("The (optional) cereal library was marked as required, but wasn't found.")
        else ()
            seqan3_config_print ("Optional dependency:        Cereal not found.")
        endif ()
    endif ()
endif ()

# ----------------------------------------------------------------------------
# ZLIB dependency
# ----------------------------------------------------------------------------

if (NOT SEQAN3_NO_ZLIB)
    find_package (ZLIB QUIET)
endif ()

if (ZLIB_FOUND)
    set (SEQAN3_LIBRARIES ${SEQAN3_LIBRARIES} ${ZLIB_LIBRARIES})
    set (SEQAN3_DEPENDENCY_INCLUDE_DIRS ${SEQAN3_DEPENDENCY_INCLUDE_DIRS} ${ZLIB_INCLUDE_DIRS})
    set (SEQAN3_DEFINITIONS ${SEQAN3_DEFINITIONS} "-DSEQAN3_HAS_ZLIB=1")
    seqan3_config_print ("Optional dependency:        ZLIB-${ZLIB_VERSION_STRING} found.")
else ()
    seqan3_config_print ("Optional dependency:        ZLIB not found.")
endif ()

# ----------------------------------------------------------------------------
# BZip2 dependency
# ----------------------------------------------------------------------------

if (NOT SEQAN3_NO_BZIP2)
    find_package (BZip2 QUIET)
endif ()

if (NOT ZLIB_FOUND AND BZIP2_FOUND)
    # NOTE (marehr): iostream_bzip2 uses the type `uInt`, which is defined by
    # `zlib`. Therefore, `bzip2` will cause a ton of errors without `zlib`.
    message (AUTHOR_WARNING "Disabling BZip2 [which was successfully found], "
                            "because ZLIB was not found. BZip2 depends on ZLIB.")
    unset (BZIP2_FOUND)
endif ()

if (BZIP2_FOUND)
    set (SEQAN3_LIBRARIES ${SEQAN3_LIBRARIES} ${BZIP2_LIBRARIES})
    set (SEQAN3_DEPENDENCY_INCLUDE_DIRS ${SEQAN3_DEPENDENCY_INCLUDE_DIRS} ${BZIP2_INCLUDE_DIRS})
    set (SEQAN3_DEFINITIONS ${SEQAN3_DEFINITIONS} "-DSEQAN3_HAS_BZIP2=1")
    seqan3_config_print ("Optional dependency:        BZip2-${BZIP2_VERSION_STRING} found.")
else ()
    seqan3_config_print ("Optional dependency:        BZip2 not found.")
endif ()

# ----------------------------------------------------------------------------
# System dependencies
# ----------------------------------------------------------------------------

# librt
if ((${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
    OR (${CMAKE_SYSTEM_NAME} STREQUAL "kFreeBSD")
    OR (${CMAKE_SYSTEM_NAME} STREQUAL "GNU"))
    set (SEQAN3_LIBRARIES ${SEQAN3_LIBRARIES} rt)
endif ()

# libexecinfo -- implicit
check_include_file_cxx (execinfo.h _SEQAN3_HAVE_EXECINFO)
mark_as_advanced (_SEQAN3_HAVE_EXECINFO)
if (_SEQAN3_HAVE_EXECINFO)
    seqan3_config_print ("Optional dependency:        libexecinfo found.")
    if ((${CMAKE_SYSTEM_NAME} STREQUAL "FreeBSD") OR (${CMAKE_SYSTEM_NAME} STREQUAL "OpenBSD"))
        set (SEQAN3_LIBRARIES ${SEQAN3_LIBRARIES} execinfo elf)
    endif ()
else ()
    seqan3_config_print ("Optional dependency:        libexecinfo not found.")
endif ()

# ----------------------------------------------------------------------------
# Perform compilability test of platform.hpp (tests some requirements)
# ----------------------------------------------------------------------------

set (CXXSTD_TEST_SOURCE "#include <seqan3/core/platform.hpp>
                         int main() {}")

# using try_compile instead of check_cxx_source_compiles to capture output in case of failure
file (WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src.cxx" "${CXXSTD_TEST_SOURCE}\n")

try_compile (SEQAN3_PLATFORM_TEST #
             ${CMAKE_BINARY_DIR}
             ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src.cxx
             CMAKE_FLAGS "-DCOMPILE_DEFINITIONS:STRING=${CMAKE_CXX_FLAGS} ${SEQAN3_CXX_FLAGS}"
                         "-DINCLUDE_DIRECTORIES:STRING=${CMAKE_INCLUDE_PATH};${SEQAN3_INCLUDE_DIR};${SEQAN3_DEPENDENCY_INCLUDE_DIRS}"
             COMPILE_DEFINITIONS ${SEQAN3_DEFINITIONS}
             LINK_LIBRARIES ${SEQAN3_LIBRARIES}
             OUTPUT_VARIABLE SEQAN3_PLATFORM_TEST_OUTPUT)

if (SEQAN3_PLATFORM_TEST)
    seqan3_config_print ("SeqAn3 platform.hpp build:  passed.")
else ()
    seqan3_config_error ("SeqAn3 platform.hpp build:  failed!\n\
                        ${SEQAN3_PLATFORM_TEST_OUTPUT}")
endif ()

# ----------------------------------------------------------------------------
# Finish find_package call
# ----------------------------------------------------------------------------

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

# propagate SEQAN3_INCLUDE_DIR into SEQAN3_INCLUDE_DIRS
set (SEQAN3_INCLUDE_DIRS ${SEQAN3_INCLUDE_DIR} ${SEQAN3_DEPENDENCY_INCLUDE_DIRS})

# ----------------------------------------------------------------------------
# Export targets
# ----------------------------------------------------------------------------

if (SEQAN3_FOUND AND NOT TARGET seqan3::seqan3)
    separate_arguments (SEQAN3_CXX_FLAGS_LIST UNIX_COMMAND "${SEQAN3_CXX_FLAGS}")

    add_library (seqan3_seqan3 INTERFACE)
    target_compile_definitions (seqan3_seqan3 INTERFACE ${SEQAN3_DEFINITIONS})
    target_compile_options (seqan3_seqan3 INTERFACE ${SEQAN3_CXX_FLAGS_LIST})
    target_link_libraries (seqan3_seqan3 INTERFACE "${SEQAN3_LIBRARIES}")
    # include seqan3/include/ as -I, because seqan3 should never produce warnings.
    target_include_directories (seqan3_seqan3 INTERFACE "${SEQAN3_INCLUDE_DIR}")
    # include everything except seqan3/include/ as -isystem, i.e.
    # a system header which suppresses warnings of external libraries.
    target_include_directories (seqan3_seqan3 SYSTEM INTERFACE "${SEQAN3_DEPENDENCY_INCLUDE_DIRS}")
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
    message ("  SEQAN3_CXX_FLAGS            ${SEQAN3_CXX_FLAGS}")
    message ("")
    message ("  SEQAN3_VERSION              ${SEQAN3_VERSION}")
    message ("  SEQAN3_VERSION_MAJOR        ${SEQAN3_VERSION_MAJOR}")
    message ("  SEQAN3_VERSION_MINOR        ${SEQAN3_VERSION_MINOR}")
    message ("  SEQAN3_VERSION_PATCH        ${SEQAN3_VERSION_PATCH}")
endif ()

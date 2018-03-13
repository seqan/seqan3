# ============================================================================
#                  SeqAn - The Library for Sequence Analysis
# ============================================================================
#
# Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
# Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of Knut Reinert or the FU Berlin nor the names of
#       its contributors may be used to endorse or promote products derived
#       from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
# OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
# DAMAGE.
# ============================================================================
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
#   C++17
#   C++ Concepts (either via Concepts TS or C++20)
#   C++ Filesystem (part of C++17 but needs extra linking on some platforms)
#   pthread
#
# SeqAn requires the following libraries:
#
#   SDSL      -- the succinct data structure library
#   Range-V3  -- Ranges Library by Eric Niebler
#
# SeqAn has the following optional dependencies:
#
#   ZLIB      -- zlib compression library
#   BZip2     -- libbz2 compression library
#   Cereal    -- Serialisation library
#   Lemon     -- Graph library
#
# If you don't wish for these to be detected (and used), you may define SEQAN3_NO_ZLIB,
# SEQAN3_NO_BZIP2, SEQAN3_NO_CEREAL and SEQAN3_NO_LEMON respectively.
#
# If you wish to require the presence of ZLIB or BZip2, just check for the module before
# finding SeqAn3, e.g. "find_package (ZLIB REQUIRED)".
# If you wish to require the presence of CEREAL, you may define SEQAN3_CEREAL.
# If you wish to require the presence of LEMON, you may define SEQAN3_LEMON.
#
# Once the search has been performed, the following variables will be set.
#
#   SEQAN3_FOUND            -- Indicate whether SeqAn was found and requirements met.
#   SeqAn3_FOUND            -- the same as SEQAN3_FOUND
#   seqan3_FOUND            -- the same as SEQAN3_FOUND
#
#   SEQAN3_VERSION_STRING   -- The version as string, e.g. "3.0.0"
#   SEQAN3_VERSION_MAJOR    -- e.g. 3
#   SEQAN3_VERSION_MINOR    -- e.g. 0
#   SEQAN3_VERSION_PATCH    -- e.g. 0
#
#   SEQAN3_INCLUDE_DIRS     -- to be passed to include_directories ()
#   SEQAN3_LIBRARIES        -- to be passed to target_link_libraries ()
#   SEQAN3_DEFINITIONS      -- to be passed to add_definitions ()
#   SEQAN3_CXX_FLAGS        -- to be added to CMAKE_CXX_FLAGS
#
# ============================================================================

cmake_minimum_required (VERSION 3.2)

# ----------------------------------------------------------------------------
# Set initial variables
# ----------------------------------------------------------------------------

set (SEQAN3_FOUND FALSE)
set (seqan3_FOUND FALSE)
set (SeqAn3_FOUND FALSE)

# work around obscure case sensitivity problems in CMake (https://cmake.org/pipermail/cmake/2009-March/027414.html)
if (DEFINED seqan3_DIR)
    set (FIND_NAME "seqan3")
elseif (DEFINED SeqAn3_DIR)
    set (FIND_NAME "SeqAn3")
elseif (DEFINED SEQAN3_DIR)
    set (FIND_NAME "SEQAN3")
else ()
    message (FATAL_ERROR "You must give \"SEQAN3\", \"SeqAn3\" or \"seqan3\" as the package name to find_package ();\n \
                         other case/combinations are not supported.")
    return ()
endif ()

# ----------------------------------------------------------------------------
# Greeter
# ----------------------------------------------------------------------------

string (ASCII 27 Esc)
set (ColourBold "${Esc}[1m")
set (ColourReset "${Esc}[m")

if (NOT ${FIND_NAME}_FIND_QUIETLY)
    message (STATUS "${ColourBold}Finding SeqAn3 and checking requirements:${ColourReset}")
endif ()

# ----------------------------------------------------------------------------
# Includes
# ----------------------------------------------------------------------------

include (FindPackageMessage)
include (CheckIncludeFileCXX)
include (CheckCXXSourceCompiles)

# ----------------------------------------------------------------------------
# Pretty printing and error handling
# ----------------------------------------------------------------------------

macro (seqan3_config_print text)
    if (NOT ${FIND_NAME}_FIND_QUIETLY)
        message (STATUS "  ${text}")
    endif ()
endmacro ()

macro (seqan3_config_error text)
    if (${FIND_NAME}_FIND_REQUIRED)
        message (FATAL_ERROR ${text})
    else ()
        if (NOT ${FIND_NAME}_FIND_QUIETLY)
            message (WARNING ${text})
        endif ()
        return ()
    endif ()
endmacro ()

# ----------------------------------------------------------------------------
# Detect if we are a clone of repository and if yes auto-add submodules
# ----------------------------------------------------------------------------

if (EXISTS "${CMAKE_CURRENT_LIST_DIR}/../build_system")
    message (STATUS "  Detected as running from a repository checkout…")

    if (IS_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/../include")
        message (STATUS "  …adding SeqAn3 include:     ${CMAKE_CURRENT_LIST_DIR}/../include")
        set (SEQAN3_INCLUDE_DIRS "${CMAKE_CURRENT_LIST_DIR}/../include" ${SEQAN3_INCLUDE_DIRS})
    endif ()

    if (EXISTS "${CMAKE_CURRENT_LIST_DIR}/../submodules")
        file (GLOB submodules ${CMAKE_CURRENT_LIST_DIR}/../submodules/*/include)
        foreach (submodule ${submodules})
            if (IS_DIRECTORY ${submodule})
                message (STATUS "  …adding submodule include:  ${submodule}")
                set (SEQAN3_INCLUDE_DIRS ${submodule} ${SEQAN3_INCLUDE_DIRS})
            endif ()
        endforeach ()
    endif ()
endif ()

# ----------------------------------------------------------------------------
# Options for CheckCXXSourceCompiles
# ----------------------------------------------------------------------------

# deactivate messages in check_*
set (CMAKE_REQUIRED_QUIET       1)
# use global variables in Check* calls
set (CMAKE_REQUIRED_INCLUDES    ${CMAKE_INCLUDE_PATH} ${SEQAN3_INCLUDE_DIRS})
set (CMAKE_REQUIRED_FLAGS       ${CMAKE_CXX_FLAGS})

# ----------------------------------------------------------------------------
# Force-deactivate optional dependencies
# ----------------------------------------------------------------------------

# Cereal is auto-detected by default, i.e. used if found, not used if not found.
# You can optionally set a hard requirement so a build fails without cereal,
# or you can force-disable cereal even if present on the system.
option (SEQAN3_CEREAL    "Require cereal and fail if not present." OFF)
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

# Lemon is auto-detected by default, i.e. used if found, not used if not found.
# You can optionally set a hard requirement so a build fails without Lemon,
# or you can force-disable Lemon even if present on the system.
option (SEQAN3_LEMON    "Require Lemon and fail if not present." OFF)
option (SEQAN3_NO_LEMON "Don't use Lemon, even if present." OFF)

if (SEQAN3_LEMON AND SEQAN3_NO_LEMON)
    # this is always a user error, therefore we always error-out, even if SeqAn is not required
    message (FATAL_ERROR "You may not specify SEQAN3_LEMON and SEQAN3_NO_LEMON at the same time.\n\
                          You can specify neither (use auto-detection), or specify either to force on/off.")
    return ()
endif ()

if (SEQAN3_LEMON)
    set (SEQAN3_DEFINITIONS ${SEQAN3_DEFINITIONS} "-DSEQAN3_WITH_LEMON=1")
elseif (SEQAN3_NO_LEMON)
    set (SEQAN3_DEFINITIONS ${SEQAN3_DEFINITIONS} "-DSEQAN3_WITH_LEMON=0")
endif ()

# These two are "opt-in", because detected by CMake
# If you want to force-require these, just do find_package (zlib REQUIRED) before find_package (seqan3)
option (SEQAN3_NO_ZLIB  "Don't use ZLIB, even if present." OFF)
option (SEQAN3_NO_BZIP2 "Don't use BZip2, even if present." OFF)

# ----------------------------------------------------------------------------
# Require C++17
# ----------------------------------------------------------------------------

set (CMAKE_REQUIRED_FLAGS_ORIGINAL ${CMAKE_REQUIRED_FLAGS})

set (CXXSTD_TEST_SOURCE
    "#if !defined (__cplusplus) || (__cplusplus < 201703L)
    #error NOCXX17
    #endif
    int main() {}")

check_cxx_source_compiles ("${CXXSTD_TEST_SOURCE}" CXX17_BUILTIN)

if (CXX17_BUILTIN)
    seqan3_config_print ("C++ Standard-17 support:    builtin")
else ()
    set (CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS_ORIGINAL} -std=c++17")

    check_cxx_source_compiles ("${CXXSTD_TEST_SOURCE}" CXX17_FLAG)

    if (CXX17_FLAG)
        seqan3_config_print ("C++ Standard-17 support:    via -std=c++17")
    else ()
        seqan3_config_error ("SeqAn3 requires C++17, but your compiler does not support it.")
    endif ()

    set (SEQAN3_CXX_FLAGS "${SEQAN3_CXX_FLAGS} -std=c++17")
endif ()

# ----------------------------------------------------------------------------
# Require C++ Concepts
# ----------------------------------------------------------------------------

set (CMAKE_REQUIRED_FLAGS_ORIGINAL ${CMAKE_REQUIRED_FLAGS})

set (CXXSTD_TEST_SOURCE
    "static_assert (__cpp_concepts >= 201507);
    int main() {}")

check_cxx_source_compiles ("${CXXSTD_TEST_SOURCE}" CONCEPTS_BUILTIN)

if (CONCEPTS_BUILTIN)
    seqan3_config_print ("C++ Concepts support:       builtin")
else ()
    set (CONCEPTS_FLAG "")

    foreach (_FLAG -std=c++20 -std=c++2a -fconcepts)
        set (CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS_ORIGINAL} ${_FLAG}")

        check_cxx_source_compiles ("${CXXSTD_TEST_SOURCE}" CONCEPTS_FLAG${_FLAG})

        if (CONCEPTS_FLAG${_FLAG})
            set (SEQAN3_CXX_FLAGS "${SEQAN3_CXX_FLAGS} ${_FLAG}")
            set (CONCEPTS_FLAG ${_FLAG})
            break ()
        endif ()
    endforeach ()

    if (CONCEPTS_FLAG)
        seqan3_config_print ("C++ Concepts support:       via ${CONCEPTS_FLAG}")
    else ()
        seqan3_config_error ("SeqAn3 requires C++ Concepts, but your compiler does not support them.")
    endif ()
endif ()

# ----------------------------------------------------------------------------
# Require C++ Filesystem
# ----------------------------------------------------------------------------

# find the correct header
check_include_file_cxx (filesystem _SEQAN3_HAVE_FILESYSTEM)
check_include_file_cxx (experimental/filesystem _SEQAN3_HAVE_EXP_FILESYSTEM)

if (_SEQAN3_HAVE_FILESYSTEM)
    seqan3_config_print ("C++ Filesystem header:      <filesystem>")

    set (CXXSTD_TEST_SOURCE
        "#include <filesystem>
        int main()
        {
            std::filesystem::path p{\"\tmp/\"};
        }")
elseif (_SEQAN3_HAVE_EXP_FILESYSTEM)
    seqan3_config_print ("C++ Filesystem header:      <experimental/filesystem>")

    set (CXXSTD_TEST_SOURCE
        "#include <experimental/filesystem>
        int main()
        {
            std::experimental::filesystem::path p{\"/tmp/\"};
        }")
else ()
    seqan3_config_error ("SeqAn3 requires C++17 filesystem support, but the filesystem header was not found.")
endif ()

# check if library is required
set (CMAKE_REQUIRED_LIBRARIES_ORIGINAL ${CMAKE_REQUIRED_LIBRARIES})

check_cxx_source_compiles ("${CXXSTD_TEST_SOURCE}" C++17FS_BUILTIN)

if (C++17FS_BUILTIN)
    seqan3_config_print ("C++ Filesystem library:     builtin")
else ()
    set (C++17FS_LIB "")

    foreach (_LIB stdc++fs)
        set (CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES_ORIGINAL} ${_LIB})

        check_cxx_source_compiles ("${CXXSTD_TEST_SOURCE}" C++17FS_LIB-l${_LIB})

        if (C++17FS_LIB-l${_LIB})
            set (SEQAN3_LIBRARIES ${SEQAN3_LIBRARIES} ${_LIB})
            set (C++17FS_LIB ${_LIB})
            break ()
        endif ()
    endforeach ()

    if (C++17FS_LIB)
        seqan3_config_print ("C++ Filesystem library:     via -l${C++17FS_LIB}")
    else ()
        seqan3_config_error ("SeqAn3 requires C++17 filesystem support, but your compiler does not offer it.")
    endif ()
endif ()

# ----------------------------------------------------------------------------
# Require Ranges and SDSL
# ----------------------------------------------------------------------------

check_include_file_cxx (range/v3/version.hpp _SEQAN3_HAVE_RANGEV3)

if (_SEQAN3_HAVE_RANGEV3)
    seqan3_config_print ("Required dependency:        Range-V3 found.")
else ()
    seqan3_config_error ("The range-v3 library is required, but wasn't found. Get it from https://github.com/ericniebler/range-v3/")
endif ()

# TODO: Doesn't have a version file, yet
# check_include_file_cxx (sdsl-lite/version.hpp _SEQAN3_HAVE_SDSL)
# if (NOT _SEQAN3_HAVE_SDSL)
#     seqan3_config_error ("The SDSL library is required, but wasn't found. Get it from https://github.com/xxsds/sdsl-lite")
# endif ()

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
# Lemon dependency is optional, but may set as required
# ----------------------------------------------------------------------------

if (NOT SEQAN3_NO_LEMON)
    check_include_file_cxx (lemon/config.h _SEQAN3_HAVE_LEMON)

    if (_SEQAN3_HAVE_LEMON)
        if (SEQAN3_LEMON)
            seqan3_config_print ("Required dependency:        Lemon found.")
        else ()
            seqan3_config_print ("Optional dependency:        Lemon found.")
        endif ()
    else ()
        if (SEQAN3_LEMON)
            seqan3_config_error ("The (optional) Lemon library was marked as required, but wasn't found.")
        else ()
            seqan3_config_print ("Optional dependency:        Lemon not found.")
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
    set (SEQAN3_LIBRARIES         ${SEQAN3_LIBRARIES}         ${ZLIB_LIBRARIES})
    set (SEQAN3_INCLUDE_DIRS      ${SEQAN3_INCLUDE_DIRS}      ${ZLIB_INCLUDE_DIRS})
    set (SEQAN3_DEFINITIONS       ${SEQAN3_DEFINITIONS}       "-DSEQAN3_HAS_ZLIB=1")
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
    set (SEQAN3_LIBRARIES         ${SEQAN3_LIBRARIES}         ${BZIP2_LIBRARIES})
    set (SEQAN3_INCLUDE_DIRS      ${SEQAN3_INCLUDE_DIRS}      ${BZIP2_INCLUDE_DIRS})
    set (SEQAN3_DEFINITIONS       ${SEQAN3_DEFINITIONS}       "-DSEQAN3_HAS_BZIP2=1")
    seqan3_config_print ("Optional dependency:        BZip2-${BZIP2_VERSION_STRING} found.")
else ()
    seqan3_config_print ("Optional dependency:        BZip2 not found.")
endif ()

# ----------------------------------------------------------------------------
# System dependencies
# ----------------------------------------------------------------------------

# librt
if ((${CMAKE_SYSTEM_NAME} STREQUAL "Linux") OR
    (${CMAKE_SYSTEM_NAME} STREQUAL "kFreeBSD") OR
    (${CMAKE_SYSTEM_NAME} STREQUAL "GNU"))
    set (SEQAN3_LIBRARIES ${SEQAN3_LIBRARIES} rt)
endif ()

# libpthread TODO (h-2): re-evaluate whether we do complete static linking in the future
if (UNIX AND NOT APPLE)
    # pthread flag (without -l) also defines macros on some platforms
    set (SEQAN3_CXX_FLAGS "${SEQAN3_CXX_FLAGS} -pthread")
    # need to explicitly link lpthread with whole-archive because Debian strips bits from the library
    # that are required for C++14 threading, resulting in broken static builds
    set (SEQAN3_LIBRARIES ${SEQAN3_LIBRARIES} "-Wl,--whole-archive -lpthread -Wl,--no-whole-archive")
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
# Find SeqAn3 include path
# ----------------------------------------------------------------------------

if (NOT SEQAN3_BASEDIR)
    find_path (SEQAN3_BASEDIR "seqan3"  PATHS ${SEQAN3_INCLUDE_DIRS})
endif ()

mark_as_advanced (SEQAN3_BASEDIR)

# find include directory
if (SEQAN3_BASEDIR)
    set (SEQAN3_INCLUDE_DIRS ${SEQAN3_BASEDIR} ${SEQAN3_INCLUDE_DIRS})
    seqan3_config_print ("SeqAn3 include dir found:   ${SEQAN3_BASEDIR}")
else ()
    seqan3_config_error ("SeqAn3 include directory could not be found.")
endif ()

# ----------------------------------------------------------------------------
# Find SeqAn3 version.hpp and extract version
# ----------------------------------------------------------------------------

set (_SEQAN3_VERSION_HPP "")
foreach (_INCLUDE_DIR ${SEQAN3_INCLUDE_DIRS})
    get_filename_component (_SEQAN3_VERSION_HPP "${_INCLUDE_DIR}/seqan3/version.hpp" ABSOLUTE)
    if (EXISTS ${_SEQAN3_VERSION_HPP})
        break ()
    endif ()
endforeach ()

set (_SEQAN3_VERSION_IDS MAJOR MINOR PATCH)

# set to 0.0.0 ny default
foreach (_ID ${_SEQAN3_VERSION_IDS})
    set (_SEQAN3_VERSION_${_ID} "0")
endforeach ()

# Exit if version.hpp is not found
if (EXISTS "${_SEQAN3_VERSION_HPP}")
    seqan3_config_print ("SeqAn3 version.hpp found:   ${_SEQAN3_VERSION_HPP}")

    foreach (_ID ${_SEQAN3_VERSION_IDS})
        file (STRINGS ${_SEQAN3_VERSION_HPP} _VERSION_${_ID} REGEX ".*SEQAN3_VERSION_${_ID}.*")
        string (REGEX REPLACE ".*SEQAN3_VERSION_${_ID}[ |\t]+([0-9a-zA-Z]+).*" "\\1" _SEQAN3_VERSION_${_ID} ${_VERSION_${_ID}})
    endforeach ()
else ()
    seqan3_config_error ("seqan3/version.hpp not found. Broken or missing install of SeqAn3.")
endif ()

set (_SEQAN3_VERSION_STRING "${_SEQAN3_VERSION_MAJOR}.${_SEQAN3_VERSION_MINOR}.${_SEQAN3_VERSION_PATCH}")

# Cache results.
set (SEQAN3_VERSION_MAJOR "${_SEQAN3_VERSION_MAJOR}" CACHE INTERNAL "SeqAn major version.")
set (SEQAN3_VERSION_MINOR "${_SEQAN3_VERSION_MINOR}" CACHE INTERNAL "SeqAn minor version.")
set (SEQAN3_VERSION_PATCH "${_SEQAN3_VERSION_PATCH}" CACHE INTERNAL "SeqAn patch version.")
set (SEQAN3_VERSION_STRING "${_SEQAN3_VERSION_STRING}" CACHE INTERNAL "SeqAn version string.")

if (SEQAN3_VERSION_STRING EQUAL 0.0.0)
    seqan3_config_error ("SeqAn version could not be detected from file. Aborting.")
elseif ((SEQAN3_VERSION_STRING VERSION_LESS 3.0.0) OR (NOT SEQAN3_VERSION_STRING VERSION_LESS 4.0.0))
    seqan3_config_error ("SeqAn version detected is ${SEQAN3_VERSION_STRING}; this is not SeqAn3!")
else ()
    seqan3_config_print ("SeqAn3 version detected:    ${SEQAN3_VERSION_STRING}")
endif ()

# ----------------------------------------------------------------------------
# Perform compilability test of platform.hpp (tests some requirements)
# ----------------------------------------------------------------------------

set (CXXSTD_TEST_SOURCE
     "#include <seqan3/core/platform.hpp>
     int main() {}")

# using try_compile instead of check_cxx_source_compiles to capture output in case of failure
file (WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src.cxx" "${CXXSTD_TEST_SOURCE}\n")

try_compile (SEQAN3_PLATFORM_TEST
             ${CMAKE_BINARY_DIR}
             ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src.cxx
             CMAKE_FLAGS         "-DCOMPILE_DEFINITIONS:STRING=${CMAKE_CXX_FLAGS} ${SEQAN3_CXX_FLAGS}"
                                 "-DINCLUDE_DIRECTORIES:STRING=${CMAKE_INCLUDE_PATH};${SEQAN3_INCLUDE_DIRS}"
             COMPILE_DEFINITIONS ${SEQAN3_DEFINITIONS}
             LINK_LIBRARIES      ${SEQAN3_LIBRARIES}
             OUTPUT_VARIABLE     SEQAN3_PLATFORM_TEST_OUTPUT)

if (SEQAN3_PLATFORM_TEST)
    seqan3_config_print ("SeqAn3 platform.hpp build:  passed.")
else ()
    seqan3_config_error ("SeqAn3 platform.hpp build:  failed!\n\
                        ${SEQAN3_PLATFORM_TEST_OUTPUT}")
endif ()

# ----------------------------------------------------------------------------
# We made it!
# ----------------------------------------------------------------------------

set (SEQAN3_FOUND TRUE)
set (seqan3_FOUND TRUE)
set (SeqAn3_FOUND TRUE)

# ----------------------------------------------------------------------------
# Print Variables
# ----------------------------------------------------------------------------

if (NOT ${FIND_NAME}_FIND_QUIETLY)
    message (STATUS "${ColourBold}Found SeqAn3:${ColourReset} ${SEQAN3_BASEDIR}/seqan3 (found version \"${SEQAN3_VERSION_STRING}\")")
endif ()

if (SEQAN3_FIND_DEBUG)
  message ("Result for ${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt")
  message ("")
  message ("  CMAKE_BUILD_TYPE            ${CMAKE_BUILD_TYPE}")
  message ("  CMAKE_SOURCE_DIR            ${CMAKE_SOURCE_DIR}")
  message ("  CMAKE_INCLUDE_PATH          ${CMAKE_INCLUDE_PATH}")
  message ("  SEQAN3_BASEDIR              ${SEQAN3_BASEDIR}")
  message ("")
  message ("  ${FIND_NAME}_FOUND                ${${FIND_NAME}_FOUND}")
  message ("  SEQAN3_HAS_ZLIB             ${ZLIB_FOUND}")
  message ("  SEQAN3_HAS_BZIP2            ${BZIP2_FOUND}")
  message ("")
  message ("  SEQAN3_INCLUDE_DIRS         ${SEQAN3_INCLUDE_DIRS}")
  message ("  SEQAN3_LIBRARIES            ${SEQAN3_LIBRARIES}")
  message ("  SEQAN3_DEFINITIONS          ${SEQAN3_DEFINITIONS}")
  message ("  SEQAN3_CXX_FLAGS            ${SEQAN3_CXX_FLAGS}")
  message ("")
  message ("  SEQAN3_VERSION_STRING       ${SEQAN3_VERSION_STRING}")
  message ("  SEQAN3_VERSION_MAJOR        ${SEQAN3_VERSION_MAJOR}")
  message ("  SEQAN3_VERSION_MINORG       ${SEQAN3_VERSION_MINOR}")
  message ("  SEQAN3_VERSION_PATCH        ${SEQAN3_VERSION_PATCH}")
endif ()

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

cmake_minimum_required (VERSION 3.2)

### Find doxygen and dependency to DOT tool
message (STATUS "Searching for doxygen.")
find_package (Doxygen REQUIRED)

if (NOT ${DOXYGEN_FOUND})
    message (FATAL_ERROR "Could not find doxygen. Not building documentation.")
endif ()

if (NOT ${DOXYGEN_DOT_FOUND})
    message (STATUS "Could not find dot tool. Disabling dot support.")
    set (SEQAN3_DOXYGEN_HAVE_DOT "NO")
else ()
    message (STATUS "Found dot tool. Enabling dot support.")
    set (SEQAN3_DOXYGEN_HAVE_DOT "YES")
endif ()

### Configure doc/developer targets.
set(seqan3_doxyfile_in ${SEQAN3_DOXYGEN_INPUT_DIR}/seqan3_doxygen_cfg.in)

option(SEQAN3_USER_DOC "Create build target and test for user documentation." ON)
option(SEQAN3_DEV_DOC "Create build target and test for developer documentation." ON)

if (SEQAN3_USER_DOC)
    message (STATUS "Configuring user doc.")

    set (SEQAN3_DOXYGEN_OUTPUT_DIR "${PROJECT_BINARY_DIR}/user_doc")
    set (SEQAN3_DOXYGEN_INCLUDE_DIR "${SEQAN3_INCLUDE_DIR}")
    set (SEQAN3_DOXYGEN_EXCLUDE_SYMBOLS "detail") #/""
    set (SEQAN3_DOXYGEN_PREDEFINED_NDEBUG "-NDEBUG") #/""
    set (SEQAN3_DOXYGEN_ENABLED_SECTIONS "") #/"DEV"
    set (SEQAN3_DOXYGEN_EXTRACT_PRIVATE "NO") #/"YES":
    set (seqan3_doxyfile_user ${PROJECT_BINARY_DIR}/seqan3_doxygen_cfg_user)

    configure_file (${seqan3_doxyfile_in} ${seqan3_doxyfile_user})

    add_custom_target(doc_user
                      COMMAND ${DOXYGEN_EXECUTABLE} ${seqan3_doxyfile_user}
                      WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
                      COMMENT "Generating user API documentation with Doxygen"
                      VERBATIM)
endif (SEQAN3_USER_DOC)

if (SEQAN3_DEV_DOC)
    message(STATUS "Configuring devel doc.")

    set(SEQAN3_DOXYGEN_OUTPUT_DIR "${PROJECT_BINARY_DIR}/devel_doc")
    set(SEQAN3_DOXYGEN_INCLUDE_DIR "${SEQAN3_INCLUDE_DIR}")
    set(SEQAN3_DOXYGEN_EXCLUDE_SYMBOLS "")
    set(SEQAN3_DOXYGEN_PREDEFINED_NDEBUG "")
    set(SEQAN3_DOXYGEN_ENABLED_SECTIONS "DEV")
    set(SEQAN3_DOXYGEN_EXTRACT_PRIVATE "YES")
    set(seqan3_doxyfile_devel ${PROJECT_BINARY_DIR}/seqan3_doxygen_cfg_devel)

    configure_file(${seqan3_doxyfile_in} ${seqan3_doxyfile_devel})

    add_custom_target(doc_devel
                      COMMAND ${DOXYGEN_EXECUTABLE} ${seqan3_doxyfile_devel}
                      WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
                      COMMENT "Generating developer API documentation with Doxygen"
                      VERBATIM)
                      message (STATUS "Add devel doc test.")
endif ()

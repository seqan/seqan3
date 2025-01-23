# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

### Find doxygen and dependency to DOT tool
message (STATUS "Searching for doxygen.")
find_package (Doxygen REQUIRED)

if (NOT ${DOXYGEN_FOUND})
    message (FATAL_ERROR "Could not find doxygen. Not building documentation.")
endif ()

if (NOT ${DOXYGEN_HAVE_DOT})
    message (STATUS "Could not find dot tool. Disabling dot support.")
    set (SEQAN3_DOXYGEN_HAVE_DOT "NO")
else ()
    message (STATUS "Found dot tool. Enabling dot support.")
    set (SEQAN3_DOXYGEN_HAVE_DOT "YES")
endif ()

### Number of threads to use for dot. Doxygen's default is 0 (all threads).
set (SEQAN3_DOXYGEN_DOT_NUM_THREADS "0")

### Configure doc/developer targets.
set (SEQAN3_DOXYGEN_SOURCE_DIR "${SEQAN3_CLONE_DIR}")
set (SEQAN3_DOXYFILE_IN ${SEQAN3_DOXYGEN_INPUT_DIR}/seqan3_doxygen_cfg.in)
set (SEQAN3_FOOTER_HTML_IN ${SEQAN3_DOXYGEN_INPUT_DIR}/seqan3_footer.html.in)
# DoxygenLayout.xml.in is created by seqan3-doxygen-layout.cmake
set (SEQAN3_LAYOUT_IN ${CMAKE_CURRENT_BINARY_DIR}/DoxygenLayout.xml.in)

option (SEQAN3_USER_DOC "Create build target and test for user documentation." ON)
option (SEQAN3_DEV_DOC "Create build target and test for developer documentation." ON)

### Download and extract cppreference-doxygen-web.tag.xml for std:: documentation links
# SEQAN3_DOXYGEN_STD_TAGFILE can be used to point to an existing tag file (cppreference-doxygen-web.tag.xml).
# If SEQAN3_DOXYGEN_STD_TAGFILE is set by the user and the file exists, it will be copied.
# If SEQAN3_DOXYGEN_STD_TAGFILE is not set by the user, or it is set by the user, but the file does not exist,
# the tag file will be downloaded.
set (SEQAN3_DEFAULT_DOXYGEN_STD_TAGFILE "${PROJECT_BINARY_DIR}/cppreference-doxygen-web.tag.xml")
set (SEQAN3_DOXYGEN_STD_TAGFILE
     "${SEQAN3_DEFAULT_DOXYGEN_STD_TAGFILE}"
     CACHE STRING "Path to cppreference-doxygen-web.tag.xml")
if (NOT EXISTS "${SEQAN3_DOXYGEN_STD_TAGFILE}" OR SEQAN3_DOXYGEN_STD_TAGFILE STREQUAL
                                                  "${SEQAN3_DEFAULT_DOXYGEN_STD_TAGFILE}")
    message (STATUS "Tag file will be fetched.")
    # Reset path in case it was set from the outside, but does not exist.
    set (SEQAN3_DOXYGEN_STD_TAGFILE "${SEQAN3_DEFAULT_DOXYGEN_STD_TAGFILE}")
    include (ExternalProject)
    # When updating, check whether warnings in SEQAN3_TEST_DOXYGEN_FAIL_ON_WARNINGS are gone when removing sed filter.
    ExternalProject_Add (
        download-cppreference-doxygen-web-tag
        URL "https://github.com/PeterFeicht/cppreference-doc/releases/download/v20241110/html-book-20241110.tar.xz"
        URL_HASH SHA256=431e80862eb70fd4793a60d7d3b6c13c8605284978f9ea0529572e8fd1562cc6
        TLS_VERIFY ON
        DOWNLOAD_DIR "${PROJECT_BINARY_DIR}"
        DOWNLOAD_NAME "html-book.tar.xz"
        DOWNLOAD_NO_EXTRACT YES
        BINARY_DIR "${PROJECT_BINARY_DIR}"
        BUILD_BYPRODUCTS "${SEQAN3_DEFAULT_DOXYGEN_STD_TAGFILE}"
        CONFIGURE_COMMAND /bin/sh -c "xzcat html-book.tar.xz | tar -xf - cppreference-doxygen-web.tag.xml"
        BUILD_COMMAND rm "html-book.tar.xz"
        INSTALL_COMMAND "")
else ()
    message (STATUS "Copying existing tag file: ${SEQAN3_DOXYGEN_STD_TAGFILE}")
    # Copy tag file such that it is present in the built documentation. This is useful if the documentation is
    # subsequently deployed to a server.
    add_custom_target (download-cppreference-doxygen-web-tag)
    add_custom_command (TARGET download-cppreference-doxygen-web-tag
                        POST_BUILD
                        COMMAND ${CMAKE_COMMAND} -E copy "${SEQAN3_DOXYGEN_STD_TAGFILE}"
                                "${SEQAN3_DEFAULT_DOXYGEN_STD_TAGFILE}"
                        BYPRODUCTS "${SEQAN3_DEFAULT_DOXYGEN_STD_TAGFILE}")
    set (SEQAN3_DOXYGEN_STD_TAGFILE "${SEQAN3_DEFAULT_DOXYGEN_STD_TAGFILE}")
endif ()

### TEST HELPER
add_test (NAME cppreference-doxygen-web-tag COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target
                                                    download-cppreference-doxygen-web-tag)

# doxygen does not show any warnings (doxygen prints warnings / errors to cerr)
# Second line filters warnings from tag file.
# Note: Because the commands are line-wise, CMake will insert a semicolon between them.
#       If this is changed to be a single line, the semicolon must be manually inserted.
set (SEQAN3_TEST_DOXYGEN_FAIL_ON_WARNINGS
     "${DOXYGEN_EXECUTABLE} -q > doxygen.cout 2> doxygen.cerr"
     "sed -i '/documented symbol '\\''T std::experimental::erase'\\'' was not declared or defined\\./d; /documented symbol '\\''T std::experimental::erase_if'\\'' was not declared or defined\\./d' \"doxygen.cerr\""
     "cat \"doxygen.cerr\""
     "test ! -s \"doxygen.cerr\""
     CACHE INTERNAL "The doxygen test command")

### install helper

# make sure that prefix path is /usr/local/share/doc/seqan3/
if (NOT DEFINED CMAKE_SIZEOF_VOID_P)
    # we need this to suppress GNUInstallDirs AUTHOR_WARNING:
    #   CMake Warning (dev) at /usr/share/cmake-3.19/Modules/GNUInstallDirs.cmake:223 (message):
    #     Unable to determine default CMAKE_INSTALL_LIBDIR directory because no
    #     target architecture is known.  Please enable at least one language before
    #     including GNUInstallDirs.
    set (CMAKE_SIZEOF_VOID_P 8)
endif ()
include (GNUInstallDirs) # this is needed to prefix the install paths

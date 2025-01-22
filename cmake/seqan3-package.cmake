# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

# This file describes how SeqAn3 will be packaged.

set (CPACK_GENERATOR "TXZ")

set (CPACK_PACKAGE_VERSION "${SEQAN3_VERSION}")
set (CPACK_PACKAGE_VENDOR "seqan")
# A description of the project, used in places such as the introduction screen of CPack-generated Windows installers.
# set (CPACK_PACKAGE_DESCRIPTION_FILE "") # TODO
set (CPACK_PACKAGE_CHECKSUM "SHA256")
set (CPACK_PACKAGE_ICON "${SEQAN3_CLONE_DIR}/test/documentation/seqan_logo.svg")
set (CPACK_RESOURCE_FILE_LICENSE "${SEQAN3_CLONE_DIR}/LICENSE.md")
set (CPACK_RESOURCE_FILE_README "${SEQAN3_CLONE_DIR}/README.md")

# Already being called on source package, i.e. CPM is already downloaded.
if (NOT CPM_DOWNLOAD_LOCATION)
    set (CPM_DOWNLOAD_LOCATION "${SEQAN3_CLONE_DIR}/cmake/CPM.cmake")
else ()
    CPMGetPackage (use_ccache)
endif ()

configure_file ("${SEQAN3_CLONE_DIR}/cmake/cpack_install.cmake.in" "${CMAKE_CURRENT_BINARY_DIR}/cpack_install.cmake"
                @ONLY)
set (CPACK_INSTALL_SCRIPT "${CMAKE_CURRENT_BINARY_DIR}/cpack_install.cmake")

# Source Package
set (CPACK_SOURCE_GENERATOR "TXZ")
# Next 2 lines would overwrite copy in cpack_install.cmake
list (APPEND CPACK_SOURCE_IGNORE_FILES "/cmake/CPM.cmake")
list (APPEND CPACK_SOURCE_IGNORE_FILES "/test/cmake/seqan3_require_ccache.cmake")
list (APPEND CPACK_SOURCE_IGNORE_FILES "/\.git($|/)")
list (APPEND CPACK_SOURCE_IGNORE_FILES "/\.github/")
list (APPEND CPACK_SOURCE_IGNORE_FILES "/\.vscode/")
list (APPEND CPACK_SOURCE_IGNORE_FILES "/build/")
list (APPEND CPACK_SOURCE_IGNORE_FILES "/submodules/")

include (CPack)

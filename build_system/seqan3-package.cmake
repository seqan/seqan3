# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
# -----------------------------------------------------------------------------------------------------

# This file describes how SeqAn3 will be packaged.

cmake_minimum_required (VERSION 3.7)

set (CPACK_GENERATOR "ZIP" "TXZ")

set (CPACK_PACKAGE_VENDOR "seqan")
# A description of the project, used in places such as the introduction screen of CPack-generated Windows installers.
# set (CPACK_PACKAGE_DESCRIPTION_FILE "") # TODO
set (CPACK_PACKAGE_CHECKSUM "SHA256")
set (CPACK_PACKAGE_ICON "${SEQAN3_CLONE_DIR}/test/documentation/seqan_logo.png")
set (CPACK_RESOURCE_FILE_LICENSE "${SEQAN3_CLONE_DIR}/LICENSE.md")
set (CPACK_RESOURCE_FILE_README "${SEQAN3_CLONE_DIR}/README.md")

include (CPack)

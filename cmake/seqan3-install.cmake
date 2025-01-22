# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

# This file describes where and which parts of SeqAn3 should be installed to.

include (GNUInstallDirs)

# install documentation files in /share/doc
install (FILES "${SEQAN3_CLONE_DIR}/CHANGELOG.md" #
               "${SEQAN3_CLONE_DIR}/CODE_OF_CONDUCT.md" #
               "${SEQAN3_CLONE_DIR}/CONTRIBUTING.md" #
               "${SEQAN3_CLONE_DIR}/LICENSE.md" #
               "${SEQAN3_CLONE_DIR}/README.md"
         TYPE DOC)

# install cmake files in /share/cmake
install (FILES "${SEQAN3_CLONE_DIR}/cmake/seqan3-config.cmake" "${SEQAN3_CLONE_DIR}/cmake/seqan3-config-version.cmake"
         DESTINATION "${CMAKE_INSTALL_DATADIR}/cmake/seqan3")

# install seqan3 header files in /include/seqan3
install (DIRECTORY "${SEQAN3_INCLUDE_DIR}/seqan3" TYPE INCLUDE)

install (DIRECTORY "${SEQAN3_SDSL_INCLUDE_DIR}/sdsl" DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/seqan3/vendor")
install (DIRECTORY "${SEQAN3_CEREAL_INCLUDE_DIR}/cereal" DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/seqan3/vendor")

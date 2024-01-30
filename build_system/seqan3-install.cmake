# SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

# This file describes where and which parts of SeqAn3 should be installed to.

cmake_minimum_required (VERSION 3.14)

include (GNUInstallDirs)

# install documentation files in /share/doc
install (FILES "${SEQAN3_CLONE_DIR}/CHANGELOG.md" #
               "${SEQAN3_CLONE_DIR}/CODE_OF_CONDUCT.md" #
               "${SEQAN3_CLONE_DIR}/CONTRIBUTING.md" #
               "${SEQAN3_CLONE_DIR}/LICENSE.md" #
               "${SEQAN3_CLONE_DIR}/README.md"
         TYPE DOC)

# install cmake files in /share/cmake
install (FILES "${SEQAN3_CLONE_DIR}/build_system/seqan3-config.cmake"
               "${SEQAN3_CLONE_DIR}/build_system/seqan3-config-version.cmake"
         DESTINATION "${CMAKE_INSTALL_DATADIR}/cmake/seqan3")

# install seqan3 header files in /include/seqan3
install (DIRECTORY "${SEQAN3_INCLUDE_DIR}/seqan3" TYPE INCLUDE)

# install submodule header files, e.g. all external dependencies in /home/user/seqan3/submodules/*,
# in /include/seqan3/submodules/<submodule>/include
foreach (submodule_dir ${SEQAN3_DEPENDENCY_INCLUDE_DIRS})
    # e.g. submodule_dir: (1) /home/user/seqan3/submodules/sdsl-lite/include or (2) /usr/include
    # strip /home/user/seqan3/submodules/ and /include part.
    file (RELATIVE_PATH submodule "${SEQAN3_SUBMODULES_DIR}/submodules" "${submodule_dir}/..")
    # submodule is either a single module name, like sdsl-lite or a relative path to a folder ../../../usr
    # skip relative folders and only keep submodules that reside in the submodules folder
    if (NOT submodule MATCHES "^\\.\\.") # skip relative folders
        install (DIRECTORY "${submodule_dir}" DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/seqan3/submodules/${submodule}")
    endif ()
endforeach ()

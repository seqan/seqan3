# SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0

# CPM Package Lock
# This file should be committed to version control

# URL/GIT_TAG may be annotated with a branch name
# This is needed for https://github.com/seqan/actions/tree/main/update_cpm_package_lock

# The first argument of CPMDeclarePackage can be freely chosen and is used as argument in CPMGetPackage.
# The NAME argument should be package name that would also be used in a find_package call.
# Ideally, both are the same, which might not always be possible: https://github.com/cpm-cmake/CPM.cmake/issues/603
# This is needed to support CPM_USE_LOCAL_PACKAGES

# cereal
set (SEQAN3_CEREAL_VERSION 1.3.2)
CPMDeclarePackage (cereal
                   NAME cereal
                   VERSION ${SEQAN3_CEREAL_VERSION}
                   GITHUB_REPOSITORY USCiLab/cereal
                   DOWNLOAD_ONLY TRUE
                   QUIET YES)
# sdsl-lite
# Use URL download of the commit archive such that we do not clone submodules
# Package name is still sdsl (name as v2 at xxsds/sdsl), but sdsl-lite is not currently being packaged
# To avoid accidentally using the older sdsl, NAME is set to sdsl-lite
set (SEQAN3_SDSL_VERSION 14cd017027ea742353fc5b500d1cb1d95896b77e)
CPMDeclarePackage (sdsl-lite
                   NAME sdsl-lite
                   URL https://github.com/xxsds/sdsl-lite/archive/${SEQAN3_SDSL_VERSION}.tar.gz # master
                   DOWNLOAD_ONLY YES
                   QUIET YES)
# benchmark
set (SEQAN3_BENCHMARK_VERSION 1.9.0)
CPMDeclarePackage (benchmark
                   NAME benchmark
                   VERSION ${SEQAN3_BENCHMARK_VERSION}
                   GITHUB_REPOSITORY google/benchmark
                   SYSTEM TRUE
                   OPTIONS "BENCHMARK_ENABLE_TESTING OFF" "BENCHMARK_ENABLE_WERROR OFF"
                           "CMAKE_MESSAGE_LOG_LEVEL WARNING")
# googletest
set (SEQAN3_GOOGLETEST_VERSION 1.15.2)
CPMDeclarePackage (googletest
                   NAME GTest
                   VERSION ${SEQAN3_GOOGLETEST_VERSION}
                   GITHUB_REPOSITORY google/googletest
                   SYSTEM TRUE
                   OPTIONS "BUILD_GMOCK OFF" "INSTALL_GTEST OFF" "CMAKE_CXX_STANDARD 20"
                           "CMAKE_MESSAGE_LOG_LEVEL WARNING")
# doxygen-awesome
set (SEQAN3_DOXYGEN_AWESOME_VERSION 2.3.4)
CPMDeclarePackage (doxygen_awesome
                   NAME doxygen_awesome
                   VERSION ${SEQAN3_DOXYGEN_AWESOME_VERSION}
                   GITHUB_REPOSITORY jothepro/doxygen-awesome-css
                   DOWNLOAD_ONLY TRUE
                   QUIET YES)
# seqan2
set (SEQAN3_SEQAN2_VERSION 1f1575a9e05c06a6b9efb50d23ce60e540cbaaf2)
CPMDeclarePackage (seqan
                   NAME seqan
                   GIT_TAG ${SEQAN3_SEQAN2_VERSION} # main
                   GITHUB_REPOSITORY seqan/seqan
                   DOWNLOAD_ONLY YES
                   QUIET YES)
# use_ccache
set (SEQAN3_USE_CCACHE_VERSION d2a54ef555b6fc2d496a4c9506dbeb7cf899ce37)
CPMDeclarePackage (use_ccache
                   NAME use_ccache
                   GIT_TAG ${SEQAN3_USE_CCACHE_VERSION} # main
                   GITHUB_REPOSITORY seqan/cmake-scripts
                   SOURCE_SUBDIR ccache
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE)

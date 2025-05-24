# SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0

# CPM Package Lock
# This file should be committed to version control

# URL/GIT_TAG may be annotated with a branch name
# This is needed for https://github.com/seqan/actions/tree/main/update_cpm_package_lock

# The first argument of CPMDeclarePackage can be freely chosen and is used as argument in CPMGetPackage.
# The NAME argument should be package name that would also be used in a find_package call.
# Ideally, both are the same, which might not always be possible: https://github.com/cpm-cmake/CPM.cmake/issues/603
# This is needed to support CPM_USE_LOCAL_PACKAGES

# Each package has a (project-prefixed) version variable, which allows changing the version of a package without
# changing the package lock file.
# This is useful for packaging, where there might only be an older version, e.g. of googletest, available.
# Note that the variable has to be a cache variable to work properly, but not a forced cache variable.
# A clean reconfigure, i.e. deleting CMakeCache.txt, is needed to update the default version in an existing
# build directory.

# cmake-format: off

# cereal
set (SEQAN3_CEREAL_VERSION 1.3.2 CACHE STRING "")
CPMDeclarePackage (cereal
                   NAME cereal
                   VERSION ${SEQAN3_CEREAL_VERSION}
                   GITHUB_REPOSITORY USCiLab/cereal
                   SYSTEM TRUE
                   OPTIONS "JUST_INSTALL_CEREAL ON" "CMAKE_MESSAGE_LOG_LEVEL WARNING")
# benchmark
set (SEQAN3_BENCHMARK_VERSION 1.9.4 CACHE STRING "")
CPMDeclarePackage (benchmark
                   NAME benchmark
                   VERSION ${SEQAN3_BENCHMARK_VERSION}
                   GITHUB_REPOSITORY google/benchmark
                   SYSTEM TRUE
                   OPTIONS "BENCHMARK_ENABLE_TESTING OFF" "BENCHMARK_ENABLE_WERROR OFF"
                           "CMAKE_MESSAGE_LOG_LEVEL WARNING")
# googletest
set (SEQAN3_GOOGLETEST_VERSION 1.17.0 CACHE STRING "")
CPMDeclarePackage (googletest
                   NAME GTest
                   VERSION ${SEQAN3_GOOGLETEST_VERSION}
                   GITHUB_REPOSITORY google/googletest
                   SYSTEM TRUE
                   OPTIONS "BUILD_GMOCK OFF" "INSTALL_GTEST OFF" "CMAKE_MESSAGE_LOG_LEVEL WARNING")
# doxygen-awesome
set (SEQAN3_DOXYGEN_AWESOME_VERSION 2.3.4 CACHE STRING "")
CPMDeclarePackage (doxygen_awesome
                   NAME doxygen_awesome
                   VERSION ${SEQAN3_DOXYGEN_AWESOME_VERSION}
                   GITHUB_REPOSITORY jothepro/doxygen-awesome-css
                   DOWNLOAD_ONLY TRUE
                   QUIET TRUE)
# seqan2
set (SEQAN3_SEQAN2_VERSION 2.5.0 CACHE STRING "")
CPMDeclarePackage (seqan
                   NAME seqan
                   VERSION ${SEQAN3_SEQAN2_VERSION}
                   GIT_TAG seqan-v${SEQAN3_SEQAN2_VERSION}
                   GITHUB_REPOSITORY seqan/seqan
                   DOWNLOAD_ONLY TRUE
                   QUIET TRUE)
# use_ccache
set (SEQAN3_USE_CCACHE_VERSION d2a54ef555b6fc2d496a4c9506dbeb7cf899ce37 CACHE STRING "")
CPMDeclarePackage (use_ccache
                   NAME use_ccache
                   GIT_TAG ${SEQAN3_USE_CCACHE_VERSION} # main
                   GITHUB_REPOSITORY seqan/cmake-scripts
                   SOURCE_SUBDIR ccache
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE)

# cmake-format: on

# SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0

# CPM Package Lock
# This file should be committed to version control

# cereal
set (SEQAN3_CEREAL_VERSION 1.3.2)
CPMDeclarePackage (cereal
                   NAME cereal
                   VERSION ${SEQAN3_CEREAL_VERSION}
                   GITHUB_REPOSITORY USCiLab/cereal
                   DOWNLOAD_ONLY TRUE)
# sdsl-lite
set (SEQAN3_SDSL_VERSION 14cd017027ea742353fc5b500d1cb1d95896b77e)
CPMDeclarePackage (sdsl-lite
                   NAME sdsl-lite
                   GIT_TAG ${SEQAN3_SDSL_VERSION}
                   GITHUB_REPOSITORY xxsds/sdsl-lite
                   DOWNLOAD_ONLY YES)
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
                   NAME googletest
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
                   DOWNLOAD_ONLY TRUE)

# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

# list all search places
# NOTE: this can be enabled globally by -DSEQAN3_EXTERNAL_PROJECT_FIND_DEBUG_MODE=1 to help debug search paths
# set (CMAKE_FIND_DEBUG_MODE ${SEQAN3_EXTERNAL_PROJECT_FIND_DEBUG_MODE})

message (STATUS "CMAKE_COMMAND: ${CMAKE_COMMAND}")
message (STATUS "CMAKE_VERSION: ${CMAKE_VERSION}")
message (STATUS "CMAKE_SYSTEM_NAME: ${CMAKE_SYSTEM_NAME}")

message (STATUS "search paths: (see https://cmake.org/cmake/help/latest/command/find_package.html#search-procedure)")
message (STATUS "1) CMAKE_FIND_USE_PACKAGE_ROOT_PATH: ${CMAKE_FIND_USE_CMAKE_PATH}")
message (STATUS "   SeqAn3_ROOT: ${SeqAn3_ROOT}")

message (STATUS "2) CMAKE_FIND_USE_CMAKE_PATH: ${CMAKE_FIND_USE_CMAKE_PATH}")
message (STATUS "   CMAKE_PREFIX_PATH: ${CMAKE_PREFIX_PATH}")
message (STATUS "   CMAKE_FRAMEWORK_PATH: ${CMAKE_FRAMEWORK_PATH}")
message (STATUS "   CMAKE_APPBUNDLE_PATH: ${CMAKE_APPBUNDLE_PATH}")

message (STATUS "3) CMAKE_FIND_USE_CMAKE_ENVIRONMENT_PATH: ${CMAKE_FIND_USE_CMAKE_ENVIRONMENT_PATH}")
message (STATUS "   SeqAn3_DIR: ${SeqAn3_DIR}")

# 4) with hints. See point 4 in https://cmake.org/cmake/help/latest/command/find_package.html#search-procedure.
# There is no output.

message (STATUS "5) CMAKE_FIND_USE_SYSTEM_ENVIRONMENT_PATH: ${CMAKE_FIND_USE_SYSTEM_ENVIRONMENT_PATH}")

message (STATUS "6) CMAKE_FIND_USE_PACKAGE_REGISTRY: ${CMAKE_FIND_USE_PACKAGE_REGISTRY}")

message (STATUS "7) CMAKE_FIND_USE_CMAKE_SYSTEM_PATH: ${CMAKE_FIND_USE_CMAKE_SYSTEM_PATH}")
message (STATUS "   CMAKE_SYSTEM_PREFIX_PATH: ${CMAKE_SYSTEM_PREFIX_PATH}")
message (STATUS "   CMAKE_SYSTEM_FRAMEWORK_PATH: ${CMAKE_SYSTEM_FRAMEWORK_PATH}")
message (STATUS "   CMAKE_SYSTEM_APPBUNDLE_PATH: ${CMAKE_SYSTEM_APPBUNDLE_PATH}")

message (STATUS "8) CMAKE_FIND_USE_SYSTEM_PACKAGE_REGISTRY: ${CMAKE_FIND_USE_SYSTEM_PACKAGE_REGISTRY}")

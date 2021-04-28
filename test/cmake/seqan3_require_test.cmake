# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
# -----------------------------------------------------------------------------------------------------

# Exposes the google-test targets `gtest` and `gtest_main`.
macro (seqan3_require_test_old gtest_git_tag)
    set (gtest_project_args ${SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS})
    list (APPEND gtest_project_args "-DBUILD_GMOCK=0")

    # force that libraries are installed to `lib/`, because GNUInstallDirs might install it into `lib64/`
    list (APPEND gtest_project_args "-DCMAKE_INSTALL_LIBDIR=${PROJECT_BINARY_DIR}/lib/")

    # google sets CMAKE_DEBUG_POSTFIX = "d"
    set (gtest_main_path "${PROJECT_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}gtest_main${CMAKE_STATIC_LIBRARY_SUFFIX}")
    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
        set (gtest_main_path "${PROJECT_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}gtest_maind${CMAKE_STATIC_LIBRARY_SUFFIX}")
    endif ()

    set (gtest_path "${PROJECT_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}gtest${CMAKE_STATIC_LIBRARY_SUFFIX}")
    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
        set (gtest_path "${PROJECT_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}gtestd${CMAKE_STATIC_LIBRARY_SUFFIX}")
    endif ()

    include (ExternalProject)
    ExternalProject_Add (
        gtest_project
        PREFIX gtest_project
        GIT_REPOSITORY "https://github.com/google/googletest.git"
        GIT_TAG "${gtest_git_tag}"
        SOURCE_DIR "${SEQAN3_TEST_CLONE_DIR}"
        CMAKE_ARGS "${gtest_project_args}"
        BUILD_BYPRODUCTS "${gtest_main_path}" "${gtest_path}"
        UPDATE_DISCONNECTED ${SEQAN3_TEST_BUILD_OFFLINE}
    )
    unset (gtest_project_args)

    add_library (gtest_main STATIC IMPORTED)
    add_dependencies (gtest_main gtest_project)
    set_target_properties (gtest_main PROPERTIES IMPORTED_LOCATION "${gtest_main_path}")

    add_library (gtest STATIC IMPORTED)
    add_dependencies (gtest gtest_main)
    set_target_properties (gtest PROPERTIES IMPORTED_LOCATION "${gtest_path}")
    set_property (TARGET gtest APPEND PROPERTY INTERFACE_LINK_LIBRARIES "pthread")

    unset(gtest_main_path)
    unset(gtest_path)
endmacro ()

macro (seqan3_require_test)
    enable_testing ()

    set (gtest_git_tag "252ce9c52d304659eff6be558209c811b7191963") # 26-04-2021

    if (NOT CMAKE_VERSION VERSION_LESS 3.14)
        message (STATUS "Fetch googletest:")

        include (FetchContent)
        FetchContent_Declare (
            gtest_fetch_content
            GIT_REPOSITORY "https://github.com/google/googletest.git"
            GIT_TAG "${gtest_git_tag}"
        )
        option (BUILD_GMOCK "" OFF)
        FetchContent_MakeAvailable(gtest_fetch_content)
    else ()
        message (STATUS "Use googletest as external project:")

        seqan3_require_test_old ("${gtest_git_tag}")
    endif ()

    add_custom_target (gtest_build DEPENDS gtest_main gtest)
endmacro ()

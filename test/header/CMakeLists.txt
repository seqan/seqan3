# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
# -----------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.7)
project (seqan3_header_test CXX)

include (../seqan3-test.cmake)

option (SEQAN3_FULL_HEADER_TEST "Test seqan3 headers as well as the headers of external libraries" OFF)

# We compile each header twice in separate compilation units, where the second compilation omits the header guard (check
# for cyclic includes). Each alone is sufficient to test that the header is functional, but both are needed to check for
# link errors, which can happen if the header accidentally defines a variable, e.g. a global or class static member.
# Furthermore this tests that header guards are working by including the same header twice.
#
# example invocation:
#     seqan3_header_test (seqan3 "<path>/include/seqan3" "test.hpp|/folder/|test2.hpp")
#
# \param component        The component name, will create the target `${component}_header_test`
# \param header_base_path The base path to the header files
# \param exclude_regex    A regular expression on the header file paths that excludes them from the header test. For
#                         regex syntax see https://cmake.org/cmake/help/v3.15/command/string.html#regular-expressions.
#
# \sa Modified version from Bio-Formats
# https://github.com/openmicroscopy/bioformats/blob/d3bb33eeda23e81f78fd25f658bfc14a4363805f/cpp/cmake/HeaderTest.cmake#L81-L113
#
# \sa https://en.cppreference.com/w/cpp/language/storage_duration#external_linkage
#
# ======== Detailed Explanation ========
#
# `test/header/generate_header_source.cmake` generates for each given header file two `.cpp` files.
#
# Let header be `include/seqan3/search/fm_index/bi_fm_index_cursor.hpp` and component be `seqan3`.
#
# * `<build_dir>/seqan3_header_test/seqan3/search/fm_index/bi_fm_index_cursor.hpp-header-guard.cpp`
#   which checks if the header implements a header guard.
#
#   This will be checked by including the header twice, if the header guard is missing the second include will produce
#   redeclaration errors, i.e.
#
#   ```c++
#   #include <seqan3/search/fm_index/bi_fm_index_cursor.hpp>
#   #include <seqan3/search/fm_index/bi_fm_index_cursor.hpp> // redeclaration errors if bi_fm_index_cursor.hpp has no header-guards
#   ```
#
# * `<build_dir>/seqan3_header_test/seqan3/search/fm_index/bi_fm_index_cursor.hpp-no-self-include.cpp`
#   which checks if the header will not be included itself (cyclic include). For this we copy the content of the header
#   file into the source file and if any included header includes the header itself again, we would get redeclaration
#   errors, i.e.
#   ```c++
#   // This is file bi_fm_index_cursor.hpp without #pragma once
#
#   // ...
#   #include <seqan3/search/fm_index/fm_index.hpp> // redeclaration errors if fm_index.hpp includes bi_fm_index_cursor.hpp
#   ```
#
# Each `.cpp` file generates an object file which will be linked into one big binary. This ensures that we don't have
# any leaking symbols (e.g. https://en.cppreference.com/w/cpp/language/storage_duration#external_linkage).
macro (seqan3_header_test component header_base_path exclude_regex)
    set (target "${component}_header_test")

    # finding all *.hpp files relative from the current directory (e.g. /test/)
    # The resulting list is normalized to `header_base_path` that means concatenating
    # "${header_base_path}/header_files[i]" will result in an absolute path to the file
    #
    # Example output:
    #   seqan3/alphabet/adaptation/all.hpp
    #   seqan3/alphabet/adaptation/char.hpp
    #   seqan3/alphabet/adaptation/concept.hpp
    #   seqan3/alphabet/adaptation/uint.hpp
    #   seqan3/alphabet/all.hpp
    #   seqan3/alphabet/dna5_detail.hpp <- will be filtered out
    #   ....
    seqan3_test_files (header_files "${header_base_path}" "*.hpp;*.h")

    # filter out headers
    if (NOT ";${exclude_regex};" STREQUAL ";;")
        list (FILTER header_files EXCLUDE REGEX "${exclude_regex}")
    endif()

    file (WRITE "${PROJECT_BINARY_DIR}/${target}.cpp" "")
    add_executable (${target} ${PROJECT_BINARY_DIR}/${target}.cpp)
    target_link_libraries (${target} seqan3::test::header)
    add_test (NAME "header/${target}" COMMAND ${target})

    foreach (header ${header_files})
        seqan3_test_component (header_test_name "${header}" TEST_NAME)
        seqan3_test_component (header_target_name "${header}" TARGET_UNIQUE_NAME)

        foreach (header_sub_test "header-guard" "no-self-include")
            set (header_target_source "${PROJECT_BINARY_DIR}/${target}_files/${header_test_name}.hpp-${header_sub_test}.cpp")
            set (header_target "${target}--${header_target_name}-${header_sub_test}")

            string (REPLACE "-" "__" header_test_name_safe "${target}, ${header_target}")

            # we use add_custom_command to detect changes to a header file, which will update the generated source file
            add_custom_command (OUTPUT "${header_target_source}"
                                COMMAND "${CMAKE_COMMAND}"
                                        "-DHEADER_FILE_ABSOLUTE=${header_base_path}/${header}"
                                        "-DHEADER_FILE_INCLUDE=${header}"
                                        "-DHEADER_TARGET_SOURCE=${header_target_source}"
                                        "-DHEADER_TEST_NAME_SAFE=${header_test_name_safe}"
                                        "-DHEADER_COMPONENT=${component}"
                                        "-DHEADER_SUB_TEST=${header_sub_test}"
                                        "-P"
                                        "${CMAKE_CURRENT_SOURCE_DIR}/generate_header_source.cmake"
                                DEPENDS "${header_base_path}/${header}"
                                        "${CMAKE_CURRENT_SOURCE_DIR}/generate_header_source.cmake")

            add_library (${header_target} OBJECT "${header_target_source}")

            # NOTE: a simple target_link_libraries (${header_target} seqan3::test::header)
            # is not possible for OBJECT libraries before cmake 3.12
            if (CMAKE_VERSION VERSION_LESS 3.12)
                target_compile_options (${header_target} PRIVATE $<TARGET_PROPERTY:seqan3::test::header,INTERFACE_COMPILE_OPTIONS>)
                target_compile_definitions (${header_target} PRIVATE $<TARGET_PROPERTY:seqan3::test::header,INTERFACE_COMPILE_DEFINITIONS>)
                target_include_directories (${header_target} PRIVATE $<TARGET_PROPERTY:seqan3::test::header,INTERFACE_INCLUDE_DIRECTORIES>)
                add_dependencies (${header_target} gtest gbenchmark)
            else ()
                target_link_libraries (${header_target} seqan3::test::header)
            endif ()

            target_sources (${target} PRIVATE $<TARGET_OBJECTS:${header_target}>)
        endforeach ()
    endforeach ()

    unset (target)
    unset (header_files)
    unset (header_test_name)
    unset (header_test_name_safe)
    unset (header_target_name)
    unset (header_target_source)
    unset (header_target)
endmacro ()

seqan3_require_ccache ()
seqan3_require_benchmark ()
seqan3_require_test ()

# note: seqan3/std/* will not be tested, because the source files don't have any file extension
# note: seqan3/version.hpp is one of the only header that is not required to have a seqan3/core/platform.hpp include
seqan3_header_test (seqan3 "${SEQAN3_CLONE_DIR}/include" "seqan3/version.hpp")
seqan3_header_test (seqan3_test "${SEQAN3_CLONE_DIR}/test/include" "")

if (SEQAN3_FULL_HEADER_TEST)

    # contains deprecated headers that produce an error when included
    seqan3_header_test (range-v3 "${SEQAN3_CLONE_DIR}/submodules/range-v3/include" "/detail/")
    target_compile_definitions (range-v3_header_test INTERFACE RANGES_DISABLE_DEPRECATED_WARNINGS=1)

    # not self-contained headers; error: extra ‘;’ [-Werror=pedantic]
    # seqan3_header_test (lemon "${SEQAN3_CLONE_DIR}/submodules/lemon/include" "")

    # not complete self-contained headers
    seqan3_header_test (cereal "${SEQAN3_CLONE_DIR}/submodules/cereal/include" "/external/|polymorphic_impl\.hpp")

    seqan3_header_test (sdsl-lite "${SEQAN3_CLONE_DIR}/submodules/sdsl-lite/include" "")

endif ()

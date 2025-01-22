# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.20...3.31)
project (seqan3_test_snippet CXX)

include (../seqan3-test.cmake)
include (../cmake/diagnostics/list_unused_snippets.cmake)

CPMGetPackage (googletest)

option (SEQAN3_GENERATE_SNIPPETS "Whether seqan3 snippets should be generated and overwritten within the source tree."
        ON)

add_library (snippet_main snippet_main.cpp)
target_link_libraries (snippet_main PUBLIC seqan3::test gtest)
target_compile_options (snippet_main PUBLIC "-Wno-unused-variable" "-Wno-unused-but-set-variable")

macro (seqan3_snippet test_name_prefix snippet snippet_base_path)
    seqan3_test_component (snippet_target_name "${snippet}" TARGET_NAME)
    seqan3_test_component (snippet_test_name "${snippet}" TEST_NAME)
    seqan3_test_component (snippet_target_path "${snippet}" TARGET_PATH)
    set (target "${snippet_target_name}_snippet")

    add_executable (${target} "${snippet_base_path}/${snippet}")
    target_link_libraries (${target} PUBLIC seqan3::test snippet_main)
    set_target_properties (${target} PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                                                "${CMAKE_CURRENT_BINARY_DIR}/${snippet_target_path}")
    collect_used_snippet (${target})

    set (snippet_test_target "${test_name_prefix}/${snippet_test_name}_snippet")
    add_test (NAME "${snippet_test_target}" COMMAND ${target})

    # disable version checker, as it interferes with comparing the snippet output
    set_tests_properties ("${snippet_test_target}" PROPERTIES ENVIRONMENT SEQAN3_NO_VERSION_CHECK=0)

    set (snippet_compare_test_target "${test_name_prefix}/${snippet_test_name}_snippet_cmp_output")

    add_test (NAME "${snippet_compare_test_target}"
              COMMAND ${CMAKE_COMMAND} -DTARGET_FILE=$<TARGET_FILE:${target}>
                      -DSOURCE_FILE=${snippet_base_path}/${snippet} #
                      -P "${CMAKE_CURRENT_SOURCE_DIR}/compare_snippet_output.cmake")

    # disable version checker, as it interferes with comparing the snippet output
    set_tests_properties ("${snippet_compare_test_target}" PROPERTIES ENVIRONMENT SEQAN3_NO_VERSION_CHECK=0)

    unset (snippet_target_name)
    unset (snippet_test_name)
    unset (snippet_target_path)
    unset (target)
endmacro ()

macro (seqan3_snippets test_name_prefix snippet_base_path)
    seqan3_test_files (snippet_files "${snippet_base_path}" "*.cpp")
    foreach (snippet ${snippet_files})
        seqan3_snippet ("${test_name_prefix}" "${snippet}" "${snippet_base_path}")
    endforeach ()
endmacro ()

seqan3_snippets ("snippet" "${CMAKE_CURRENT_SOURCE_DIR}")
seqan3_snippets ("doc/snippet" "${CMAKE_CURRENT_SOURCE_DIR}/../../doc")

# These are the tests that are known to have non-deterministic output.
set_tests_properties (
    "snippet/alignment/pairwise/parallel_align_pairwise_with_callback_snippet_cmp_output"
    "snippet/io/views/async_input_buffer_snippet_cmp_output"
    "snippet/release/3.0.2_user-callback-in-alignment-and-search_snippet_cmp_output"
    "snippet/utility/container/aligned_allocator_snippet_cmp_output"
    PROPERTIES SKIP_RETURN_CODE 1)

list_unused_snippets ("${CMAKE_CURRENT_SOURCE_DIR}")
list_unused_snippets ("${CMAKE_CURRENT_SOURCE_DIR}/../../doc")

if (SEQAN3_GENERATE_SNIPPETS)
    include (seqan3_generate_snippet)

    foreach (source_snippet IN
             ITEMS "test/snippet/alphabet/nucleotide/@target_alphabet@_implicit_conversion_from_@source_alphabet@.cpp.in"
                   "test/snippet/alphabet/nucleotide/@target_alphabet@_implicit_conversion_from_@source_alphabet@_inherit.cpp.in"
                   "test/snippet/alphabet/nucleotide/@target_alphabet@_implicit_conversion_from_@source_alphabet@_vector.cpp.in"
                   "test/snippet/alphabet/nucleotide/@target_alphabet@_implicit_conversion_from_@source_alphabet@_views.cpp.in"
                   "doc/fragments/@target_alphabet@_implicit_conversion_from_@source_alphabet@.md.in")
        seqan3_generate_snippet ("${source_snippet}" -Dtarget_alphabet=dna4 -Dsource_alphabet=rna4)
        seqan3_generate_snippet ("${source_snippet}" -Dtarget_alphabet=rna4 -Dsource_alphabet=dna4)
        seqan3_generate_snippet ("${source_snippet}" -Dtarget_alphabet=dna5 -Dsource_alphabet=rna5)
        seqan3_generate_snippet ("${source_snippet}" -Dtarget_alphabet=rna5 -Dsource_alphabet=dna5)
        seqan3_generate_snippet ("${source_snippet}" -Dtarget_alphabet=dna15 -Dsource_alphabet=rna15)
        seqan3_generate_snippet ("${source_snippet}" -Dtarget_alphabet=rna15 -Dsource_alphabet=dna15)
    endforeach ()

    foreach (target_alphabet IN
             ITEMS dna4
                   rna4
                   dna5
                   rna5
                   dna15
                   rna15
                   dna16sam
                   dna3bs)
        set (source_snippet "test/snippet/alphabet/nucleotide/@target_alphabet@_char_literal.cpp.in")
        seqan3_generate_snippet ("${source_snippet}" -Dtarget_alphabet=${target_alphabet})

        set (source_snippet "test/snippet/alphabet/nucleotide/@target_alphabet@_literal.cpp.in")
        seqan3_generate_snippet ("${source_snippet}" -Dtarget_alphabet=${target_alphabet})
    endforeach ()
endif ()

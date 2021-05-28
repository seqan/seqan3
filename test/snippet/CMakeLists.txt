# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
# -----------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.7)
project (seqan3_test_snippet CXX)

include (../seqan3-test.cmake)

macro (seqan3_snippet test_name_prefix snippet snippet_base_path)
    seqan3_test_component (snippet_target_name "${snippet}" TARGET_NAME)
    seqan3_test_component (snippet_test_name "${snippet}" TEST_NAME)
    seqan3_test_component (snippet_target_path "${snippet}" TARGET_PATH)
    set (target "${snippet_target_name}_snippet")

    add_executable (${target} "${snippet_base_path}/${snippet}")
    target_link_libraries (${target} seqan3::test::unit)
    target_compile_options (${target} PUBLIC "-Wno-unused")
    set_target_properties(${target}
        PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${snippet_target_path}"
    )

    add_test (NAME "${test_name_prefix}/${snippet_test_name}_snippet" COMMAND ${target})

    unset (snippet_target_name)
    unset (snippet_test_name)
    unset (snippet_target_path)
    unset (target)
endmacro ()

macro (seqan3_snippets test_name_prefix snippet_base_path)
    seqan3_test_files (snippet_files "${snippet_base_path}" "*.cpp")
    foreach (snippet ${snippet_files})
        seqan3_snippet("${test_name_prefix}" "${snippet}" "${snippet_base_path}")
    endforeach ()
endmacro ()

seqan3_require_ccache ()
seqan3_require_test ()

seqan3_snippets ("doc/snippet" "${CMAKE_SOURCE_DIR}/../../doc")
seqan3_snippets ("snippet" "${CMAKE_SOURCE_DIR}")

include (seqan3_generate_snippet)

foreach (source_snippet IN ITEMS "test/snippet/alphabet/nucleotide/@target_alphabet@_implicit_conversion_from_@source_alphabet@.cpp.in"
                                 "test/snippet/alphabet/nucleotide/@target_alphabet@_implicit_conversion_from_@source_alphabet@_inherit.cpp.in"
                                 "test/snippet/alphabet/nucleotide/@target_alphabet@_implicit_conversion_from_@source_alphabet@_vector.cpp.in"
                                 "test/snippet/alphabet/nucleotide/@target_alphabet@_implicit_conversion_from_@source_alphabet@_views.cpp.in"
                                 "doc/fragments/@target_alphabet@_implicit_conversion_from_@source_alphabet@.hpp.in")
    seqan3_generate_snippet("${source_snippet}" -Dtarget_alphabet=dna4 -Dsource_alphabet=rna4)
    seqan3_generate_snippet("${source_snippet}" -Dtarget_alphabet=rna4 -Dsource_alphabet=dna4)
    seqan3_generate_snippet("${source_snippet}" -Dtarget_alphabet=dna5 -Dsource_alphabet=rna5)
    seqan3_generate_snippet("${source_snippet}" -Dtarget_alphabet=rna5 -Dsource_alphabet=dna5)
    seqan3_generate_snippet("${source_snippet}" -Dtarget_alphabet=dna15 -Dsource_alphabet=rna15)
    seqan3_generate_snippet("${source_snippet}" -Dtarget_alphabet=rna15 -Dsource_alphabet=dna15)
endforeach ()

foreach (target_alphabet IN ITEMS dna4 rna4 dna5 rna5 dna15 rna15 dna16sam dna3bs)
    set (source_snippet "test/snippet/alphabet/nucleotide/@target_alphabet@_char_literal.cpp.in")
    seqan3_generate_snippet("${source_snippet}" -Dtarget_alphabet=${target_alphabet})

    set (source_snippet "test/snippet/alphabet/nucleotide/@target_alphabet@_literal.cpp.in")
    seqan3_generate_snippet("${source_snippet}" -Dtarget_alphabet=${target_alphabet})
endforeach ()

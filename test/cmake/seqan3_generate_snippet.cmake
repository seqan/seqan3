# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

# Generate snippets from a source_snippet
#
# Example:
# seqan3_generate_snippet("<..>/@target_alphabet@_implicit_conversion_from_@source_alphabet@.cpp.in"
#                         -Dtarget_alphabet=dna4
#                         -Dsource_alphabet=rna4)
#
# Will generate a snippet in the same directory as the source_snippet.
#
# The generated snippet
# * will be <..>/dna4_implicit_conversion_from_rna4.cpp
# * the @placeholder@ will be replaced by the specified `-Dplaceholder=value` value
# * the .in ending will be stripped
#
# In the source_snippet
# * all ${placeholder} will be replaced by the specified `-Dplaceholder=value` value
# * see `cmake`s [configure_file](https://cmake.org/cmake/help/latest/command/configure_file.html) for more detail
function (seqan3_generate_snippet source_snippet)
    if (NOT SEQAN3_CLONE_DIR)
        message (AUTHOR_WARNING "seqan3_generate_snippet can't be used if "
                                "SEQAN3_CLONE_DIR (i.e. no git checkout) is not defined.")
        return ()
    endif ()
    # e.g. source_snippet: <...>/@target_alphabet@_implicit_conversion_from_@source_alphabet@.cpp.in

    foreach (definition IN LISTS ARGN)
        # e.g. definition -Dsource_alphabet=rna4
        if (definition MATCHES "-D(.+)=(.+)")
            # e.g. CMAKE_MATCH_1: source_alphabet and CMAKE_MATCH_2: rna4
            # i.e. set ("source_alphabet" "rna4")
            set ("${CMAKE_MATCH_1}" "${CMAKE_MATCH_2}")
        endif ()
    endforeach ()

    # remove .in file ending
    # e.g. <...>/@target_alphabet@_implicit_conversion_from_@source_alphabet@.cpp.in
    # to:  <...>/@target_alphabet@_implicit_conversion_from_@source_alphabet@.cpp
    string (REGEX REPLACE ".in$" "" target_snippet "${source_snippet}")

    # substitute @variable@s
    #
    # substitute
    #   <...>/@target_alphabet@_implicit_conversion_from_@source_alphabet@.cpp
    # by e.g. target_alphabet: dna4 and source_alphabet: rna4
    #   <...>/dna4_implicit_conversion_from_rna4.cpp
    string (CONFIGURE "${target_snippet}" target_snippet @ONLY)

    # substitute source_snippet with all definitions
    configure_file ("${SEQAN3_CLONE_DIR}/${source_snippet}" "${SEQAN3_CLONE_DIR}/${target_snippet}")
endfunction ()

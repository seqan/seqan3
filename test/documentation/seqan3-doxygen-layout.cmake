# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
# -----------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.10)

include(${SEQAN3_INCLUDE_DIR}/../test/cmake/seqan3_test_files.cmake)

# Replaces documentation entries in variable `DOXYGEN_LAYOUT`
#
# Important: the variable of the parent scope must be named `DOXYGEN_LAYOUT` for this function to work properly.
#
# Example:
# replace_in_doxygen_layout ("${SEQAN3_INCLUDE_DIR}/../doc/about/" "About")
#
# (1) Find all 'index.md' files in ${SEQAN3_INCLUDE_DIR}/../doc/about/
# (2) For each 'index.md'
#     (a) Use the very first line, e.g. `# FOOOMOO {#fooooo_bar_ref}`,
#         for the title, e.g. "FOOOMOO" and doxygen reference, e.g. "fooooo_bar_ref".
#     (b) Generate doxygen html layout entry, e.g.
#         "      <tab type=\"user\" visible=\"yes\" title=\"FOOOMOO" url=\"\\\\ref fooooo_bar_ref\" intro=\"\"/>\n"
#         and append it to a list
# (3) Replace the doxygen html layout entry list from (2) in the ${DOXYGEN_LAYOUT} input variable
#
function (replace_in_doxygen_layout doc_path doxygen_layout_tag)
    set(DOXYGEN_LAYOUT_TAG_LINE "<tab type=\"usergroup\" visible=\"yes\" title=\"${doxygen_layout_tag}\" intro=\"\">\n")
    set(DOXYGEN_LAYOUT_DOC_PAGES ${DOXYGEN_LAYOUT_TAG_LINE}) # append header line

    # iterate over all index.md
    seqan3_test_files (doc_how_to_filenames "${doc_path}" "index.md")
    foreach (doc_how_to_filename ${doc_how_to_filenames})
        set(doc_howto_filepath "${doc_path}/${doc_how_to_filename}")
        execute_process(COMMAND head -n 1 ${doc_howto_filepath} OUTPUT_VARIABLE DOC_HEADER_LINE)
        string(REGEX MATCH "^# \(.*\) {#\(.*\)}" DUMMY ${DOC_HEADER_LINE})
        set(doc_title ${CMAKE_MATCH_1})
        set(doc_ref_name ${CMAKE_MATCH_2})
        string(APPEND DOXYGEN_LAYOUT_DOC_PAGES "      <tab type=\"user\" visible=\"yes\" title=\"${doc_title}\" url=\"\\\\ref ${doc_ref_name}\" intro=\"\"/>\n")

        unset(doc_howto_filepath)
        unset(doc_title)
        unset(doc_ref_name)
    endforeach ()

    # Replace header line and appended list of doc entries with header line
    string(REGEX REPLACE "${DOXYGEN_LAYOUT_TAG_LINE}" "${DOXYGEN_LAYOUT_DOC_PAGES}" NEW_DOXYGEN_LAYOUT ${DOXYGEN_LAYOUT})

    set(DOXYGEN_LAYOUT ${NEW_DOXYGEN_LAYOUT} PARENT_SCOPE) # replace new Doxygen layout

    unset(DOXYGEN_LAYOUT_TAG_LINE)
    unset(DOXYGEN_LAYOUT_DOC_PAGES)
endfunction ()

### Add all documentation pages to DoxygenLayout.xml
### ------------------------------------------------
# Note: variable name DOXYGEN_LAYOUT must not be changed because it is directly used within `replace_in_doxygen_layout`
file(READ "${CMAKE_SOURCE_DIR}/DoxygenLayout.xml" DOXYGEN_LAYOUT)

replace_in_doxygen_layout("${SEQAN3_INCLUDE_DIR}/../doc/about/" "About")
replace_in_doxygen_layout("${SEQAN3_INCLUDE_DIR}/../doc/setup/" "Setup")
replace_in_doxygen_layout("${SEQAN3_INCLUDE_DIR}/../doc/tutorial/" "Tutorial")
replace_in_doxygen_layout("${SEQAN3_INCLUDE_DIR}/../doc/howto/" "How-To")

file(WRITE "${CMAKE_BINARY_DIR}/DoxygenLayout.xml" ${DOXYGEN_LAYOUT})

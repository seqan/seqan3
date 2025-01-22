# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

include (${SEQAN3_INCLUDE_DIR}/../test/cmake/seqan3_test_files.cmake)

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
# Optionally, `HIDE_FROM_USER` can be used, to hide certain pages from the user documentation.
# For example,
# replace_in_doxygen_layout ("${SEQAN3_INCLUDE_DIR}/../doc/setup/" "Setup" HIDE_FROM_USER "setup_tests")
# Will hide the page containing {#setup_tests} in the first line from the user documentation.
# HIDE_FROM_USER takes multiple values, e.g., `HIDE_FROM_USER "setup_tests" "setup"`
#
function (replace_in_doxygen_layout doc_path doxygen_layout_tag)
    # Parse extra arguments: <prefix> <options> <single_value> <multi_value>
    # <options> <single_value> are not used
    cmake_parse_arguments (ARGS "" "" "HIDE_FROM_USER" ${ARGN})
    # Will parse a call like
    # replace_in_doxygen_layout ("${SEQAN3_INCLUDE_DIR}/../doc/setup/" "Setup" HIDE_FROM_USER "setup_tests")
    # and put arguments of HIDE_FROM_USER in the variable ARGS_HIDE_FROM_USER

    set (DOXYGEN_LAYOUT_TAG_LINE
         "<tab type=\"usergroup\" visible=\"yes\" title=\"${doxygen_layout_tag}\" intro=\"\">\n")
    set (DOXYGEN_LAYOUT_DOC_PAGES ${DOXYGEN_LAYOUT_TAG_LINE}) # append header line

    # iterate over all index.md
    seqan3_test_files (doc_how_to_filenames "${doc_path}" "index.md")
    foreach (doc_how_to_filename ${doc_how_to_filenames})
        set (doc_howto_filepath "${doc_path}/${doc_how_to_filename}")
        execute_process (COMMAND head -n 1 ${doc_howto_filepath} OUTPUT_VARIABLE DOC_HEADER_LINE)
        string (REGEX MATCH "^# \(.*\) {#\(.*\)}" DUMMY ${DOC_HEADER_LINE})
        if ("${DUMMY}" STREQUAL "")
            message (FATAL_ERROR "Could not parse header line of ${doc_howto_filepath}.")
        endif ()
        set (doc_title ${CMAKE_MATCH_1})
        set (doc_ref_name ${CMAKE_MATCH_2})
        set (visibility "yes")
        if ("${doc_ref_name}" IN_LIST ARGS_HIDE_FROM_USER)
            set (visibility "\${SEQAN3_SHOW_DEV_DOCS}")
        endif ()
        string (APPEND
                DOXYGEN_LAYOUT_DOC_PAGES
                "      <tab type=\"user\" visible=\"${visibility}\" title=\"${doc_title}\" url=\"\\\\ref ${doc_ref_name}\" intro=\"\"/>\n"
        )

        unset (doc_howto_filepath)
        unset (doc_title)
        unset (doc_ref_name)
        unset (visibility)
    endforeach ()

    # Replace header line and appended list of doc entries with header line
    string (REGEX REPLACE "${DOXYGEN_LAYOUT_TAG_LINE}" "${DOXYGEN_LAYOUT_DOC_PAGES}" NEW_DOXYGEN_LAYOUT
                          ${DOXYGEN_LAYOUT})

    set (DOXYGEN_LAYOUT
         ${NEW_DOXYGEN_LAYOUT}
         PARENT_SCOPE) # replace new Doxygen layout

    unset (DOXYGEN_LAYOUT_TAG_LINE)
    unset (DOXYGEN_LAYOUT_DOC_PAGES)
endfunction ()

### Add all documentation pages to DoxygenLayout.xml
### ------------------------------------------------
# Note: variable name DOXYGEN_LAYOUT must not be changed because it is directly used within `replace_in_doxygen_layout`
file (READ "${SEQAN3_DOXYGEN_INPUT_DIR}/DoxygenLayout.xml" DOXYGEN_LAYOUT)

replace_in_doxygen_layout ("${SEQAN3_INCLUDE_DIR}/../doc/about/" "About")
replace_in_doxygen_layout ("${SEQAN3_INCLUDE_DIR}/../doc/setup/" "Setup" HIDE_FROM_USER "setup_tests")
replace_in_doxygen_layout ("${SEQAN3_INCLUDE_DIR}/../doc/tutorial/" "Tutorial")
replace_in_doxygen_layout ("${SEQAN3_INCLUDE_DIR}/../doc/howto/" "How-To")

file (WRITE "${CMAKE_CURRENT_BINARY_DIR}/DoxygenLayout.xml.in" ${DOXYGEN_LAYOUT})

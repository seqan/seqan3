# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

message (STATUS "Configuring user doc.")

set (SEQAN3_DOXYGEN_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}")
set (SEQAN3_DOXYGEN_EXCLUDE_SYMBOLS "seqan3::detail seqan3::simd")
set (SEQAN3_DOXYGEN_EXCLUDE_PATTERNS "*/utility/parallel/*.hpp */utility/simd/*.hpp */alignment/pairwise/policy/*.hpp")
set (SEQAN3_DOXYGEN_PREDEFINED_NDEBUG "-NDEBUG") #/""
set (SEQAN3_DOXYGEN_ENABLED_SECTIONS "") #/"DEV"
set (SEQAN3_SHOW_DEV_DOCS "NO") #/"YES"

configure_file (${SEQAN3_DOXYFILE_IN} ${SEQAN3_DOXYGEN_OUTPUT_DIR}/Doxyfile)
configure_file (${SEQAN3_FOOTER_HTML_IN} ${SEQAN3_DOXYGEN_OUTPUT_DIR}/seqan3_footer.html)
configure_file (${SEQAN3_LAYOUT_IN} ${SEQAN3_DOXYGEN_OUTPUT_DIR}/DoxygenLayout.xml)

add_custom_target (doc_usr ALL
                   COMMAND ${DOXYGEN_EXECUTABLE}
                   WORKING_DIRECTORY ${SEQAN3_DOXYGEN_OUTPUT_DIR}
                   BYPRODUCTS html/
                   DEPENDS download-cppreference-doxygen-web-tag
                   COMMENT "Generating user API documentation with Doxygen"
                   VERBATIM)

# Install doc_usr documentation in ./export folder
# make DESTDIR=export install
install (DIRECTORY "${SEQAN3_DOXYGEN_OUTPUT_DIR}/html"
         DESTINATION "${CMAKE_INSTALL_DOCDIR}"
         COMPONENT doc
         PATTERN "*.md5" EXCLUDE
         PATTERN "*.map" EXCLUDE
         PATTERN "formula.repository" EXCLUDE)

### Enable testing

enable_testing ()

message (STATUS "Add user doc test.")
add_test (NAME doc_usr_no_warnings_test COMMAND bash -c "${SEQAN3_TEST_DOXYGEN_FAIL_ON_WARNINGS}")
set_tests_properties (doc_usr_no_warnings_test PROPERTIES DEPENDS cppreference-doxygen-web-tag)

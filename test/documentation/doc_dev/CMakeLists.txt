# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

message (STATUS "Configuring devel doc.")

set (SEQAN3_DOXYGEN_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}")
set (SEQAN3_DOXYGEN_EXCLUDE_SYMBOLS "")
set (SEQAN3_DOXYGEN_EXCLUDE_PATTERNS "")
set (SEQAN3_DOXYGEN_PREDEFINED_NDEBUG "")
set (SEQAN3_DOXYGEN_ENABLED_SECTIONS "DEV")
set (SEQAN3_SHOW_DEV_DOCS "YES")

configure_file (${SEQAN3_DOXYFILE_IN} ${SEQAN3_DOXYGEN_OUTPUT_DIR}/Doxyfile)
configure_file (${SEQAN3_FOOTER_HTML_IN} ${SEQAN3_DOXYGEN_OUTPUT_DIR}/seqan3_footer.html)
configure_file (${SEQAN3_LAYOUT_IN} ${SEQAN3_DOXYGEN_OUTPUT_DIR}/DoxygenLayout.xml)

add_custom_target (doc_dev ALL
                   COMMAND ${DOXYGEN_EXECUTABLE}
                   WORKING_DIRECTORY ${SEQAN3_DOXYGEN_OUTPUT_DIR}
                   BYPRODUCTS html/
                   DEPENDS download-cppreference-doxygen-web-tag
                   COMMENT "Generating developer API documentation with Doxygen"
                   VERBATIM)

# Install doc_dev documentation in ./install_doc_dev folder
# cmake --install . --prefix install_doc_dev --component doc_dev
install (DIRECTORY "${SEQAN3_DOXYGEN_OUTPUT_DIR}/html"
         DESTINATION "${CMAKE_INSTALL_DOCDIR}"
         COMPONENT doc-dev
         EXCLUDE_FROM_ALL
         PATTERN "*.md5" EXCLUDE
         PATTERN "*.map" EXCLUDE
         PATTERN "formula.repository" EXCLUDE)

### Enable testing

enable_testing ()

message (STATUS "Add dev doc test.")
add_test (NAME doc_dev_no_warnings_test
          COMMAND bash -c "${SEQAN3_TEST_DOXYGEN_FAIL_ON_WARNINGS}"
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties (doc_dev_no_warnings_test PROPERTIES DEPENDS cppreference-doxygen-web-tag)

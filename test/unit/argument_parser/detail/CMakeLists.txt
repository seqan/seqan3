find_path (SEQAN3_TEST_LICENSE_DIR NAMES LICENSE.md HINTS "${SEQAN3_CLONE_DIR}")

add_definitions(-DSEQAN3_TEST_LICENSE_DIR="${SEQAN3_TEST_LICENSE_DIR}")
seqan3_test(format_help_test.cpp)
seqan3_test(format_html_test.cpp)
seqan3_test(format_man_test.cpp)
seqan3_test(version_check_debug_test.cpp)
seqan3_test(version_check_release_test.cpp)

# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0

seqan3_test (format_bam_test.cpp CYCLIC_DEPENDING_INCLUDES include-seqan3-io-sam_file-format_sam.hpp)
seqan3_test (format_sam_test.cpp CYCLIC_DEPENDING_INCLUDES include-seqan3-io-sam_file-format_bam.hpp)
seqan3_test (sam_file_input_test.cpp)
seqan3_test (sam_file_output_include_test.cpp)
seqan3_test (sam_file_output_test.cpp)
seqan3_test (sam_file_record_test.cpp)
seqan3_test (sam_file_seek_test.cpp)
target_compile_definitions (sam_file_seek_test PUBLIC CURRENT_SOURCE_DIR="${CMAKE_CURRENT_LIST_DIR}")
seqan3_test (sam_tag_dictionary_test.cpp)

# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0

seqan3_test (sequence_file_input_test.cpp)
seqan3_test (sequence_file_integration_test.cpp)
seqan3_test (sequence_file_integration_no_performance_test.cpp)
seqan3_test (sequence_file_output_test.cpp)
seqan3_test (sequence_file_format_embl_test.cpp)
seqan3_test (sequence_file_format_fasta_test.cpp)
seqan3_test (sequence_file_format_fasta_no_performance_test.cpp)
seqan3_test (sequence_file_format_fastq_test.cpp)
seqan3_test (sequence_file_format_fastq_no_performance_test.cpp)
seqan3_test (sequence_file_format_genbank_test.cpp)
seqan3_test (sequence_file_format_sam_test.cpp)
seqan3_test (sequence_file_record_test.cpp)
seqan3_test (sequence_file_seek_test.cpp)
target_compile_definitions (sequence_file_seek_test PUBLIC CURRENT_SOURCE_DIR="${CMAKE_CURRENT_LIST_DIR}")

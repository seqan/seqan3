# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0

add_subdirectories ()

seqan3_benchmark (format_fasta_benchmark.cpp)
seqan3_benchmark (format_fasta_no_performance_benchmark.cpp)
seqan3_benchmark (format_fastq_benchmark.cpp)
seqan3_benchmark (format_fastq_no_performance_benchmark.cpp)
seqan3_benchmark (format_sam_benchmark.cpp)
seqan3_benchmark (format_vienna_benchmark.cpp)
seqan3_benchmark (lowlevel_stream_input_benchmark.cpp)
seqan3_benchmark (lowlevel_stream_output_benchmark.cpp)
seqan3_benchmark (stream_input_benchmark.cpp)
seqan3_benchmark (stream_output_benchmark.cpp)

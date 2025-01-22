# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0

add_subdirectories ()

seqan3_benchmark (container_assignment_benchmark.cpp)
seqan3_benchmark (container_push_back_benchmark.cpp)
seqan3_benchmark (container_seq_read_benchmark.cpp)
seqan3_benchmark (container_seq_write_benchmark.cpp)
seqan3_benchmark (gap_decorator_rand_read_benchmark.cpp)
seqan3_benchmark (gap_decorator_rand_write_benchmark.cpp)
seqan3_benchmark (gap_decorator_seq_read_benchmark.cpp)
seqan3_benchmark (gap_decorator_seq_write_benchmark.cpp)

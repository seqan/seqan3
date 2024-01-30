// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/* This file includes the same performance tests as format_fasta_benchmark.cpp
 * but sets the `SEQAN3_WORKAROUND_VIEW_PERFORMANCE` variable first.
 * This assures that different definitions of functions are used.
 */
#define SEQAN3_WORKAROUND_VIEW_PERFORMANCE 0
#include "format_fasta_benchmark.cpp"

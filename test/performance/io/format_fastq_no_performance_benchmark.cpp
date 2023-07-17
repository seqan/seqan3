// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/* This file includes the same performance tests as format_fastq_benchmark.cpp
 * but sets the `SEQAN3_WORKAROUND_VIEW_PERFORMANCE` variable first.
 * This assures that different definitions of functions are used.
 */

#define SEQAN3_WORKAROUND_VIEW_PERFORMANCE 0
#include "format_fastq_benchmark.cpp"

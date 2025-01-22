// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/* This file includes the same unit tests as sequence_file_integration_test.cpp
 * but sets the `SEQAN3_WORKAROUND_VIEW_PERFORMANCE` variable first.
 * This assures that different definitions of functions are used.
 */

#define SEQAN3_WORKAROUND_VIEW_PERFORMANCE 0
#include "sequence_file_integration_test.cpp"

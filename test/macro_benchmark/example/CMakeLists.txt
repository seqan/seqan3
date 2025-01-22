# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0

seqan3_macro_benchmark (
    memcpy_macro_benchmark.cpp
    REPETITIONS 10 # optional, default: 3
    META_DATA # optional
              foo # user-specified data
              bar # user-specified data
              # seqan3_version 3.0.1 # automatically set
              # seqan3_git_version a1829272d7bafea48580af6eaac2bfbb678c3861 # automatically set
              # sdsl_version 3.0.0 # automatically set
              # sdsl_git_version 9a0d5676fd09fb8b52af214eca2d5809c9a32dbe # automatically set
)

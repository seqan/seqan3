// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/io/sequence_file/output.hpp>

int main()
{
    seqan3::sequence_file_output fout{std::cout, seqan3::format_fasta{}};
    //                          ^ no need to specify the template arguments
}

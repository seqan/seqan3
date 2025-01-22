// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <sstream>
#include <string>
#include <tuple>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/io/sam_file/output.hpp>

int main()
{
    seqan3::sam_file_output fout{std::ostringstream{}, seqan3::format_sam{}};

    auto stream_it = fout.begin();

    std::string id;
    seqan3::dna5_vector seq;

    // ...

    // assign to file iterator
    *stream_it = std::tie(seq, id);
    // is the same as:
    fout.push_back(std::tie(seq, id));
}

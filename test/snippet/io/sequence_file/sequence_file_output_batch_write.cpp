// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/io/sequence_file/output.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::sequence_file_output fout{std::ostringstream{}, seqan3::format_fasta{}};

    std::vector<std::tuple<seqan3::dna5_vector, std::string>> range{{"ACGT"_dna5, "First"},
                                                                    {"NATA"_dna5, "2nd"},
                                                                    {"GATA"_dna5, "Third"}}; // a range of "records"

    fout = range;
    // the same as:
    range | fout;
}

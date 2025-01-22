// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include <seqan3/alphabet/nucleotide/rna5.hpp>
#include <seqan3/alphabet/structure/wuss.hpp>
#include <seqan3/io/structure_file/output.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::structure_file_output fout{std::ostringstream{}, seqan3::format_vienna{}};

    std::vector<std::tuple<seqan3::rna5_vector, std::string, std::vector<seqan3::wuss51>>> range{
        {"ACGT"_rna5, "First", "...."_wuss51},
        {"NATA"_rna5, "2nd", "...."_wuss51},
        {"GATA"_rna5, "Third", "...."_wuss51}}; // a range of "records"

    range | fout;
    // the same as:
    fout = range;
}

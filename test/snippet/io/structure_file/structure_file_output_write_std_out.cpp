// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <iostream>

#include <seqan3/alphabet/nucleotide/rna4.hpp>
#include <seqan3/alphabet/structure/wuss.hpp>
#include <seqan3/io/structure_file/output.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::structure_file_output fout{std::cout, seqan3::format_vienna{}};
    //                           ^ no need to specify the template arguments

    fout.emplace_back("AACGUU"_rna4, "example_id", ".(())."_wuss51); // default order for vienna: SEQ, ID, STRUCTURE
}

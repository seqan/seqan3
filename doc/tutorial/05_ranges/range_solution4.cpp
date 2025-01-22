// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <vector>

#include <seqan3/alphabet/container/bitpacked_sequence.hpp> // include bitpacked sequence
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/argument_parser/all.hpp> // include argument parser
#include <seqan3/core/debug_stream.hpp>   // for debug_stream

int main(int argc, char ** argv)
{
    using namespace seqan3::literals;

    seqan3::argument_parser myparser("Vector-implementations-comparison", argc, argv);
    size_t size{};
    bool use_bitvector{};
    myparser.add_positional_option(size, "Size of vector");
    myparser.add_flag(use_bitvector, 'b', "bitvector", "Use bitvector instead of vector");

    try
    {
        myparser.parse();
    }
    catch (seqan3::argument_parser_error const & ext) // catch user errors
    {
        seqan3::debug_stream << "[Error] " << ext.what() << "\n";
        return -1;
    }

    if (use_bitvector)
    {
        seqan3::bitpacked_sequence<seqan3::dna4> vector;
        vector.resize(size, 'A'_dna4);
        seqan3::debug_stream << "Allocated seqan3::bitpacked_sequence<seqan3::dna4> of size " << vector.size() << '\n';
    }
    else
    {
        std::vector<seqan3::dna4> vector{size};
        seqan3::debug_stream << "Allocated std::vector<seqan3::dna4> of size " << vector.size() << '\n';
    }
}

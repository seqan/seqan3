// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

// generated from test/snippet/alphabet/nucleotide/@target_alphabet@_implicit_conversion_from_@source_alphabet@_inherit.cpp.in

//![main]
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/nucleotide/rna5.hpp>

struct my_rna5 : public seqan3::rna5
{
    // using seqan3::rna5::rna5; // uncomment to import implicit conversion shown by letter1
};

struct my_dna5 : public seqan3::dna5
{};

int main()
{
    using namespace seqan3::literals;

    // my_rna5 letter1 = 'C'_dna5; // NO automatic implicit conversion!
    // seqan3::rna5 letter2 = my_dna5{}; // seqan3::rna5 only allows implicit conversion from seqan3::dna5!
}
//![main]

#include <seqan3/utility/concept.hpp>

static_assert(seqan3::implicitly_convertible_to<seqan3::dna5, seqan3::rna5>);
static_assert(!seqan3::implicitly_convertible_to<seqan3::dna5, my_rna5>);
static_assert(!seqan3::implicitly_convertible_to<my_dna5, seqan3::rna5>);

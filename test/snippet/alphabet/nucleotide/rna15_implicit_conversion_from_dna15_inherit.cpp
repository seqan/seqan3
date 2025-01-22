// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

// generated from test/snippet/alphabet/nucleotide/@target_alphabet@_implicit_conversion_from_@source_alphabet@_inherit.cpp.in

//![main]
#include <seqan3/alphabet/nucleotide/dna15.hpp>
#include <seqan3/alphabet/nucleotide/rna15.hpp>

struct my_rna15 : public seqan3::rna15
{
    // using seqan3::rna15::rna15; // uncomment to import implicit conversion shown by letter1
};

struct my_dna15 : public seqan3::dna15
{};

int main()
{
    using namespace seqan3::literals;

    // my_rna15 letter1 = 'C'_dna15; // NO automatic implicit conversion!
    // seqan3::rna15 letter2 = my_dna15{}; // seqan3::rna15 only allows implicit conversion from seqan3::dna15!
}
//![main]

#include <seqan3/utility/concept.hpp>

static_assert(seqan3::implicitly_convertible_to<seqan3::dna15, seqan3::rna15>);
static_assert(!seqan3::implicitly_convertible_to<seqan3::dna15, my_rna15>);
static_assert(!seqan3::implicitly_convertible_to<my_dna15, seqan3::rna15>);

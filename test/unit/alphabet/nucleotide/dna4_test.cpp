// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2018, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <range/v3/view/zip.hpp>

#include "../alphabet_test_template.hpp"
#include "../alphabet_constexpr_test_template.hpp"
#include "nucleotide_test_template.hpp"

INSTANTIATE_TYPED_TEST_CASE_P(dna4, alphabet, dna4);
INSTANTIATE_TYPED_TEST_CASE_P(dna4, alphabet_constexpr, dna4);
INSTANTIATE_TYPED_TEST_CASE_P(dna4, nucleotide, dna4);

TEST(dna4, assign_char)
{
    using t = dna4;
    std::vector<char> chars
    {
        'A', 'C', 'G', 'T', 'U', 'N',
        'a', 'c', 'g', 't', 'u', 'n',
        'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V',
        'r', 'y', 's', 'w', 'k', 'm', 'b', 'd', 'h', 'v',
        '!'
    };

    std::vector<dna4> alphabets
    {
        t::A, t::C, t::G, t::T, t::U, t::A,
        t::A, t::C, t::G, t::T, t::U, t::A,
        t::A, t::C, t::C, t::A, t::G, t::A, t::C, t::A, t::A, t::A,
        t::A, t::C, t::C, t::A, t::G, t::A, t::C, t::A, t::A, t::A,
        t::A
    };

    for (auto [ chr, alp ] : ranges::view::zip(chars, alphabets))
        EXPECT_EQ((assign_char(dna4{}, chr)), alp);
}

TEST(dna4, to_char)
{
    using t = dna4;
    std::vector<char> chars
    {
        'A', 'C', 'G', 'T', 'T',
        'A'
    };

    std::vector<dna4> alphabets
    {
        t::A, t::C, t::G, t::T, t::U,
        t::UNKNOWN
    };

    for (auto [ chr, alp ] : ranges::view::zip(chars, alphabets))
        EXPECT_EQ(to_char(alp), chr);
}

TEST(dna4, literals)
{
    using namespace seqan3::literal;

    dna4_vector v;
    v.resize(5, dna4::A);
    EXPECT_EQ(v, "AAAAA"_dna4);

    std::vector<dna4> w{dna4::A, dna4::C, dna4::G, dna4::T, dna4::U, dna4::UNKNOWN};
    EXPECT_EQ(w, "ACGTTA"_dna4);
}

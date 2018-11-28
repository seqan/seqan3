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

INSTANTIATE_TYPED_TEST_CASE_P(rna4, alphabet, rna4);
INSTANTIATE_TYPED_TEST_CASE_P(rna4, alphabet_constexpr, rna4);
INSTANTIATE_TYPED_TEST_CASE_P(rna4, nucleotide, rna4);

TEST(rna4, assign_char)
{
    using t = rna4;
    std::vector<char> chars
    {
        'A', 'C', 'G', 'T', 'U', 'N',
        'a', 'c', 'g', 't', 'u', 'n',
        'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V',
        'r', 'y', 's', 'w', 'k', 'm', 'b', 'd', 'h', 'v',
        '!'
    };

    std::vector<rna4> alphabets
    {
        t::A, t::C, t::G, t::T, t::U, t::A,
        t::A, t::C, t::G, t::T, t::U, t::A,
        t::A, t::C, t::C, t::A, t::G, t::A, t::C, t::A, t::A, t::A,
        t::A, t::C, t::C, t::A, t::G, t::A, t::C, t::A, t::A, t::A,
        t::A
    };

    for (auto [ chr, alp ] : ranges::view::zip(chars, alphabets))
        EXPECT_EQ((assign_char(rna4{}, chr)), alp);
}

TEST(rna4, to_char)
{
    using t = rna4;
    std::vector<char> chars
    {
        'A', 'C', 'G', 'U', 'U',
        'A'
    };

    std::vector<rna4> alphabets
    {
        t::A, t::C, t::G, t::T, t::U,
        t::UNKNOWN
    };

    for (auto [ chr, alp ] : ranges::view::zip(chars, alphabets))
        EXPECT_EQ(to_char(alp), chr);
}

TEST(rna4, literals)
{
    using namespace seqan3::literal;

    rna4_vector v;
    v.resize(5, rna4::A);
    EXPECT_EQ(v, "AAAAA"_rna4);

    std::vector<rna4> w{rna4::A, rna4::C, rna4::G, rna4::T, rna4::U, rna4::UNKNOWN};
    EXPECT_EQ(w, "ACGUUA"_rna4);
}

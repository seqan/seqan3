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

INSTANTIATE_TYPED_TEST_CASE_P(rna15, alphabet, rna15);
INSTANTIATE_TYPED_TEST_CASE_P(rna15, alphabet_constexpr, rna15);
INSTANTIATE_TYPED_TEST_CASE_P(rna15, nucleotide, rna15);

TEST(rna15, assign_char)
{
    using t = rna15;
    std::vector<char> chars
    {
        'A', 'C', 'G', 'T', 'U', 'N',
        'a', 'c', 'g', 't', 'u', 'n',
        'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V',
        'r', 'y', 's', 'w', 'k', 'm', 'b', 'd', 'h', 'v',
        '!'
    };

    std::vector<rna15> alphabets
    {
        t::A, t::C, t::G, t::T, t::U, t::N,
        t::A, t::C, t::G, t::T, t::U, t::N,
        t::R, t::Y, t::S, t::W, t::K, t::M, t::B, t::D, t::H, t::V,
        t::R, t::Y, t::S, t::W, t::K, t::M, t::B, t::D, t::H, t::V,
        t::N
    };

    for (auto [ chr, alp ] : ranges::view::zip(chars, alphabets))
        EXPECT_EQ((assign_char(rna15{}, chr)), alp);
}

TEST(rna15, to_char)
{
    using t = rna15;
    std::vector<char> chars
    {
        'A', 'C', 'G', 'U', 'U',
        'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V',
        'N'
    };

    std::vector<rna15> alphabets
    {
        t::A, t::C, t::G, t::T, t::U,
        t::R, t::Y, t::S, t::W, t::K, t::M, t::B, t::D, t::H, t::V,
        t::UNKNOWN
    };

    for (auto [ chr, alp ] : ranges::view::zip(chars, alphabets))
        EXPECT_EQ(to_char(alp), chr);
}

TEST(rna15, literals)
{
    using namespace seqan3::literal;

    rna15_vector v;
    v.resize(5, rna15::A);
    EXPECT_EQ(v, "AAAAA"_rna15);

    std::vector<rna15> w{rna15::A, rna15::C, rna15::G, rna15::T, rna15::U, rna15::N, rna15::UNKNOWN};
    EXPECT_EQ(w, "ACGTTNN"_rna15);
}

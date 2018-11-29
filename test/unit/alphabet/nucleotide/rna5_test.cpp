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

INSTANTIATE_TYPED_TEST_CASE_P(rna5, alphabet, rna5);
INSTANTIATE_TYPED_TEST_CASE_P(rna5, alphabet_constexpr, rna5);
INSTANTIATE_TYPED_TEST_CASE_P(rna5, nucleotide, rna5);

TEST(rna5, assign_char)
{
    using t = rna5;
    std::vector<char> chars
    {
        'A', 'C', 'G', 'T', 'U', 'N',
        'a', 'c', 'g', 't', 'u', 'n',
        'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V',
        'r', 'y', 's', 'w', 'k', 'm', 'b', 'd', 'h', 'v',
        '!'
    };

    std::vector<rna5> alphabets
    {
        t::A, t::C, t::G, t::T, t::U, t::N,
        t::A, t::C, t::G, t::T, t::U, t::N,
        t::N, t::N, t::N, t::N, t::N, t::N, t::N, t::N, t::N, t::N,
        t::N, t::N, t::N, t::N, t::N, t::N, t::N, t::N, t::N, t::N,
        t::N
    };

    for (auto [ chr, alp ] : ranges::view::zip(chars, alphabets))
        EXPECT_EQ((assign_char(rna5{}, chr)), alp);
}

TEST(rna5, to_char)
{
    using t = rna5;
    std::vector<char> chars
    {
        'A', 'C', 'G', 'U', 'U',
        'N'
    };

    std::vector<rna5> alphabets
    {
        t::A, t::C, t::G, t::T, t::U,
        t::UNKNOWN
    };

    for (auto [ chr, alp ] : ranges::view::zip(chars, alphabets))
        EXPECT_EQ(to_char(alp), chr);
}

TEST(rna5, literals)
{
    using namespace seqan3::literal;

    rna5_vector v;
    v.resize(5, rna5::A);
    EXPECT_EQ(v, "AAAAA"_rna5);

    std::vector<rna5> w{rna5::A, rna5::C, rna5::G, rna5::T, rna5::U, rna5::N, rna5::UNKNOWN};
    EXPECT_EQ(w, "ACGUUNN"_rna5);
}

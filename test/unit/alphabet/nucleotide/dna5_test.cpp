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

INSTANTIATE_TYPED_TEST_CASE_P(dna5, alphabet, dna5);
INSTANTIATE_TYPED_TEST_CASE_P(dna5, alphabet_constexpr, dna5);
INSTANTIATE_TYPED_TEST_CASE_P(dna5, nucleotide, dna5);

TEST(dna5, assign_char)
{
    using t = dna5;
    std::vector<char> chars
    {
        'A', 'C', 'G', 'T', 'U', 'N',
        'a', 'c', 'g', 't', 'u', 'n',
        'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V',
        'r', 'y', 's', 'w', 'k', 'm', 'b', 'd', 'h', 'v',
        '!'
    };

    std::vector<dna5> alphabets
    {
        t::A, t::C, t::G, t::T, t::U, t::N,
        t::A, t::C, t::G, t::T, t::U, t::N,
        t::N, t::N, t::N, t::N, t::N, t::N, t::N, t::N, t::N, t::N,
        t::N, t::N, t::N, t::N, t::N, t::N, t::N, t::N, t::N, t::N,
        t::N
    };

    for (auto [ chr, alp ] : ranges::view::zip(chars, alphabets))
        EXPECT_EQ((assign_char(dna5{}, chr)), alp);
}

TEST(dna5, to_char)
{
    using t = dna5;
    std::vector<char> chars
    {
        'A', 'C', 'G', 'T', 'T',
        'N'
    };

    std::vector<dna5> alphabets
    {
        t::A, t::C, t::G, t::T, t::U,
        t::UNKNOWN
    };

    for (auto [ chr, alp ] : ranges::view::zip(chars, alphabets))
        EXPECT_EQ(to_char(alp), chr);
}

TEST(dna5, literals)
{
    using namespace seqan3::literal;

    dna5_vector v;
    v.resize(5, dna5::A);
    EXPECT_EQ(v, "AAAAA"_dna5);

    std::vector<dna5> w{dna5::A, dna5::C, dna5::G, dna5::T, dna5::U, dna5::N, dna5::UNKNOWN};
    EXPECT_EQ(w, "ACGTTNN"_dna5);
}

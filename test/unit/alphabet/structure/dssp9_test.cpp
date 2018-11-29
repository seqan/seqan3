// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2018, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <range/v3/view/zip.hpp>

#include <seqan3/alphabet/structure/dssp9.hpp>

#include "../alphabet_test_template.hpp"
#include "../alphabet_constexpr_test_template.hpp"

INSTANTIATE_TYPED_TEST_CASE_P(dssp9, alphabet, dssp9);
INSTANTIATE_TYPED_TEST_CASE_P(dssp9, alphabet_constexpr, dssp9);

// assign_char functions
TEST(dssp9, assign_char)
{
    using t = dssp9;
    std::vector<char> input
    {
        '.', '(', ')',
        ':', ',', '-', '_', '~', ';',
        '<', '>', '[', ']', '{', '}',
        'H', 'B', 'E', 'G', 'I', 'T', 'S'
    };

    std::vector<dssp9> cmp
    {
        t::X, t::X, t::X,
        t::X, t::X, t::X, t::X, t::X, t::X, t::X, t::X, t::X, t::X, t::X, t::X,
        t::H, t::B, t::E, t::G, t::I, t::T, t::S
    };

    for (auto [ ch, cm ] : ranges::view::zip(input, cmp))
        EXPECT_EQ((assign_char(dssp9{}, ch)), cm);
}

// to_char functions
TEST(dssp9, to_char)
{
    EXPECT_EQ(to_char(dssp9::H), 'H');
    EXPECT_EQ(to_char(dssp9::B), 'B');
    EXPECT_EQ(to_char(dssp9::E), 'E');
    EXPECT_EQ(to_char(dssp9::G), 'G');
    EXPECT_EQ(to_char(dssp9::I), 'I');
    EXPECT_EQ(to_char(dssp9::T), 'T');
    EXPECT_EQ(to_char(dssp9::S), 'S');
    EXPECT_EQ(to_char(dssp9::C), 'C');
    EXPECT_EQ(to_char(dssp9::X), 'X');
}

TEST(dssp9, literals)
{
    using namespace seqan3::literal;

    std::vector<dssp9> vec1;
    vec1.resize(5, dssp9::H);
    EXPECT_EQ(vec1, "HHHHH"_dssp9);

    std::vector<dssp9> vec2{dssp9::E, dssp9::H, dssp9::H, dssp9::H, dssp9::T, dssp9::G};
    EXPECT_EQ(vec2, "EHHHTG"_dssp9);
}

// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <seqan3/alphabet/structure/dssp9.hpp>
#include <seqan3/core/detail/debug_stream_range.hpp>
#include <seqan3/utility/views/zip.hpp>

#include "../alphabet_constexpr_test_template.hpp"
#include "../alphabet_test_template.hpp"
#include "../semi_alphabet_constexpr_test_template.hpp"
#include "../semi_alphabet_test_template.hpp"

using seqan3::operator""_dssp9;

INSTANTIATE_TYPED_TEST_SUITE_P(dssp9, alphabet, seqan3::dssp9, );
INSTANTIATE_TYPED_TEST_SUITE_P(dssp9, semi_alphabet_test, seqan3::dssp9, );
INSTANTIATE_TYPED_TEST_SUITE_P(dssp9, alphabet_constexpr, seqan3::dssp9, );
INSTANTIATE_TYPED_TEST_SUITE_P(dssp9, semi_alphabet_constexpr, seqan3::dssp9, );

// assign_char functions
TEST(dssp9, assign_char)
{
    std::vector<char> input
    {
        '.', '(', ')',
        ':', ',', '-', '_', '~', ';',
        '<', '>', '[', ']', '{', '}',
        'H', 'B', 'E', 'G', 'I', 'T', 'S'
    };

    std::vector<seqan3::dssp9> cmp
    {
        'X'_dssp9, 'X'_dssp9, 'X'_dssp9,
        'X'_dssp9, 'X'_dssp9, 'X'_dssp9, 'X'_dssp9, 'X'_dssp9, 'X'_dssp9,
        'X'_dssp9, 'X'_dssp9, 'X'_dssp9, 'X'_dssp9, 'X'_dssp9, 'X'_dssp9,
        'H'_dssp9, 'B'_dssp9, 'E'_dssp9, 'G'_dssp9, 'I'_dssp9, 'T'_dssp9, 'S'_dssp9
    };

    for (auto [ ch, cm ] : seqan3::views::zip(input, cmp))
        EXPECT_EQ((seqan3::assign_char_to(ch, seqan3::dssp9{})), cm);
}

// to_char functions
TEST(dssp9, to_char)
{
    EXPECT_EQ(seqan3::to_char('H'_dssp9), 'H');
    EXPECT_EQ(seqan3::to_char('B'_dssp9), 'B');
    EXPECT_EQ(seqan3::to_char('E'_dssp9), 'E');
    EXPECT_EQ(seqan3::to_char('G'_dssp9), 'G');
    EXPECT_EQ(seqan3::to_char('I'_dssp9), 'I');
    EXPECT_EQ(seqan3::to_char('T'_dssp9), 'T');
    EXPECT_EQ(seqan3::to_char('S'_dssp9), 'S');
    EXPECT_EQ(seqan3::to_char('C'_dssp9), 'C');
    EXPECT_EQ(seqan3::to_char('X'_dssp9), 'X');
}

TEST(dssp9, literals)
{

    std::vector<seqan3::dssp9> vec1;
    vec1.resize(5, 'H'_dssp9);
    EXPECT_EQ(vec1, "HHHHH"_dssp9);

    std::vector<seqan3::dssp9> vec2{'E'_dssp9, 'H'_dssp9, 'H'_dssp9, 'H'_dssp9, 'T'_dssp9, 'G'_dssp9};
    EXPECT_EQ(vec2, "EHHHTG"_dssp9);
}

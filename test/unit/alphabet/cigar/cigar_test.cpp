// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alphabet/cigar/cigar.hpp>

#include "../semi_alphabet_test_template.hpp"
#include "../semi_alphabet_constexpr_test_template.hpp"

using seqan3::get;

INSTANTIATE_TYPED_TEST_SUITE_P(cigar, semi_alphabet_test, seqan3::cigar, );
INSTANTIATE_TYPED_TEST_SUITE_P(cigar, semi_alphabet_constexpr, seqan3::cigar, );

TEST(cigar_operation, char_literal)
{
    using seqan3::operator""_cigar_operation;

    EXPECT_EQ(seqan3::to_char('M'_cigar_operation), 'M');
    EXPECT_EQ(seqan3::to_char('D'_cigar_operation), 'D');
    EXPECT_EQ(seqan3::to_char('I'_cigar_operation), 'I');
    EXPECT_EQ(seqan3::to_char('S'_cigar_operation), 'S');
    EXPECT_EQ(seqan3::to_char('H'_cigar_operation), 'H');
    EXPECT_EQ(seqan3::to_char('N'_cigar_operation), 'N');
    EXPECT_EQ(seqan3::to_char('P'_cigar_operation), 'P');
    EXPECT_EQ(seqan3::to_char('X'_cigar_operation), 'X');
    EXPECT_EQ(seqan3::to_char('='_cigar_operation), '=');
}

TEST(cigar, brace_init)
{
    using seqan3::operator""_cigar_operation;

    seqan3::cigar c1{uint32_t{223}, 'M'_cigar_operation};
    EXPECT_EQ(c1.to_string(), seqan3::small_string<11>{"223M"});
}

TEST(cigar, to_string)
{
    seqan3::cigar c1{};

    seqan3::assign_rank_to(uint32_t{223}, get<0>(c1));
    seqan3::assign_char_to('M', get<1>(c1));
    EXPECT_EQ(c1.to_string(), seqan3::small_string<11>{"223M"});
}

TEST(cigar, assign_string)
{
    seqan3::cigar c1{};
    c1.assign_string("223M");
    EXPECT_EQ(uint32_t{223}, seqan3::to_rank(get<0>(c1)));
    EXPECT_EQ('M',           get<1>(c1).to_char());
}

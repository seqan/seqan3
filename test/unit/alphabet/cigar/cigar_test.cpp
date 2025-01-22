// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/alphabet/cigar/cigar.hpp>

#include "../semi_alphabet_constexpr_test_template.hpp"
#include "../semi_alphabet_test_template.hpp"

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
    // assign from char array
    seqan3::cigar c1{};
    c1.assign_string("223M");
    EXPECT_EQ(uint32_t{223}, seqan3::to_rank(get<0>(c1)));
    EXPECT_EQ('M', get<1>(c1).to_char());

    // assign from string
    std::string s{"4S"};
    c1.assign_string(s);
    EXPECT_EQ(uint32_t{4}, seqan3::to_rank(get<0>(c1)));
    EXPECT_EQ('S', get<1>(c1).to_char());

    // assign from std::string_view
    std::string v{s};
    c1.assign_string(v);
    EXPECT_EQ(uint32_t{4}, seqan3::to_rank(get<0>(c1)));
    EXPECT_EQ('S', get<1>(c1).to_char());

    // assign from small_string
    seqan3::small_string<11> ss{"1234D"};
    c1.assign_string(ss);
    EXPECT_EQ(uint32_t{1234}, seqan3::to_rank(get<0>(c1)));
    EXPECT_EQ('D', get<1>(c1).to_char());
}

TEST(cigar, constexpr_char_literal)
{
    using seqan3::operator""_cigar_operation;
    constexpr seqan3::cigar::operation op = 'D'_cigar_operation;
    EXPECT_EQ(op.to_rank(), 1u);
}

// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/core/debug_stream/range.hpp>

#include "../alphabet_constexpr_test_template.hpp"
#include "../alphabet_test_template.hpp"
#include "../semi_alphabet_constexpr_test_template.hpp"
#include "../semi_alphabet_test_template.hpp"
#include "phred_test_template.hpp"

using seqan3::operator""_phred42;

INSTANTIATE_TYPED_TEST_SUITE_P(phred42, alphabet, seqan3::phred42, );
INSTANTIATE_TYPED_TEST_SUITE_P(phred42, semi_alphabet_test, seqan3::phred42, );
INSTANTIATE_TYPED_TEST_SUITE_P(phred42, alphabet_constexpr, seqan3::phred42, );
INSTANTIATE_TYPED_TEST_SUITE_P(phred42, semi_alphabet_constexpr, seqan3::phred42, );
INSTANTIATE_TYPED_TEST_SUITE_P(phred42, phred, seqan3::phred42, );

TEST(phred42, char_literal)
{
    EXPECT_EQ(seqan3::to_char('!'_phred42), '!');
    EXPECT_EQ(seqan3::to_char('"'_phred42), '"');
    EXPECT_EQ(seqan3::to_char('#'_phred42), '#');
    EXPECT_EQ(seqan3::to_char('$'_phred42), '$');
    EXPECT_EQ(seqan3::to_char('%'_phred42), '%');
    EXPECT_EQ(seqan3::to_char('&'_phred42), '&');
    EXPECT_EQ(seqan3::to_char('\''_phred42), '\'');
    EXPECT_EQ(seqan3::to_char('('_phred42), '(');
    EXPECT_EQ(seqan3::to_char(')'_phred42), ')');
    EXPECT_EQ(seqan3::to_char('*'_phred42), '*');
    EXPECT_EQ(seqan3::to_char('+'_phred42), '+');
    EXPECT_EQ(seqan3::to_char(','_phred42), ',');
    EXPECT_EQ(seqan3::to_char('-'_phred42), '-');
    EXPECT_EQ(seqan3::to_char('.'_phred42), '.');
    EXPECT_EQ(seqan3::to_char('/'_phred42), '/');
    EXPECT_EQ(seqan3::to_char('0'_phred42), '0');
    EXPECT_EQ(seqan3::to_char('1'_phred42), '1');
    EXPECT_EQ(seqan3::to_char('2'_phred42), '2');
    EXPECT_EQ(seqan3::to_char('3'_phred42), '3');
    EXPECT_EQ(seqan3::to_char('4'_phred42), '4');
    EXPECT_EQ(seqan3::to_char('5'_phred42), '5');
    EXPECT_EQ(seqan3::to_char('6'_phred42), '6');
    EXPECT_EQ(seqan3::to_char('7'_phred42), '7');
    EXPECT_EQ(seqan3::to_char('8'_phred42), '8');
    EXPECT_EQ(seqan3::to_char('9'_phred42), '9');
    EXPECT_EQ(seqan3::to_char(':'_phred42), ':');
    EXPECT_EQ(seqan3::to_char(';'_phred42), ';');
    EXPECT_EQ(seqan3::to_char('<'_phred42), '<');
    EXPECT_EQ(seqan3::to_char('='_phred42), '=');
    EXPECT_EQ(seqan3::to_char('>'_phred42), '>');
    EXPECT_EQ(seqan3::to_char('?'_phred42), '?');
    EXPECT_EQ(seqan3::to_char('@'_phred42), '@');
    EXPECT_EQ(seqan3::to_char('A'_phred42), 'A');
    EXPECT_EQ(seqan3::to_char('B'_phred42), 'B');
    EXPECT_EQ(seqan3::to_char('C'_phred42), 'C');
    EXPECT_EQ(seqan3::to_char('D'_phred42), 'D');
    EXPECT_EQ(seqan3::to_char('E'_phred42), 'E');
    EXPECT_EQ(seqan3::to_char('F'_phred42), 'F');
    EXPECT_EQ(seqan3::to_char('G'_phred42), 'G');
    EXPECT_EQ(seqan3::to_char('H'_phred42), 'H');
    EXPECT_EQ(seqan3::to_char('I'_phred42), 'I');
    EXPECT_EQ(seqan3::to_char('J'_phred42), 'J');
}

TEST(phred42, string_literal)
{
    std::vector<seqan3::phred42> v;
    v.resize(5, '#'_phred42);
    EXPECT_EQ(v, "#####"_phred42);

    std::vector<seqan3::phred42> w{'#'_phred42, '#'_phred42, '!'_phred42, '!'_phred42, '!'_phred42, '#'_phred42};
    EXPECT_EQ(w, "##!!!#"_phred42);
}

// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alphabet/quality/phred63.hpp>

#include "../alphabet_constexpr_test_template.hpp"
#include "../alphabet_test_template.hpp"
#include "../semi_alphabet_constexpr_test_template.hpp"
#include "../semi_alphabet_test_template.hpp"
#include "phred_test_template.hpp"

using namespace seqan3;

INSTANTIATE_TYPED_TEST_SUITE_P(phred63, alphabet_, phred63, );
INSTANTIATE_TYPED_TEST_SUITE_P(phred63, semi_alphabet_test, phred63, );
INSTANTIATE_TYPED_TEST_SUITE_P(phred63, alphabet_constexpr, phred63, );
INSTANTIATE_TYPED_TEST_SUITE_P(phred63, semi_alphabet_constexpr, phred63, );
INSTANTIATE_TYPED_TEST_SUITE_P(phred63, phred, phred63, );

TEST(phred63, char_literal)
{
    EXPECT_EQ(to_char('!'_phred63), '!');
    EXPECT_EQ(to_char('"'_phred63), '"');
    EXPECT_EQ(to_char('#'_phred63), '#');
    EXPECT_EQ(to_char('$'_phred63), '$');
    EXPECT_EQ(to_char('%'_phred63), '%');
    EXPECT_EQ(to_char('&'_phred63), '&');
    EXPECT_EQ(to_char('\''_phred63), '\'');
    EXPECT_EQ(to_char('('_phred63), '(');
    EXPECT_EQ(to_char(')'_phred63), ')');
    EXPECT_EQ(to_char('*'_phred63), '*');
    EXPECT_EQ(to_char('+'_phred63), '+');
    EXPECT_EQ(to_char(','_phred63), ',');
    EXPECT_EQ(to_char('-'_phred63), '-');
    EXPECT_EQ(to_char('.'_phred63), '.');
    EXPECT_EQ(to_char('/'_phred63), '/');
    EXPECT_EQ(to_char('0'_phred63), '0');
    EXPECT_EQ(to_char('1'_phred63), '1');
    EXPECT_EQ(to_char('2'_phred63), '2');
    EXPECT_EQ(to_char('3'_phred63), '3');
    EXPECT_EQ(to_char('4'_phred63), '4');
    EXPECT_EQ(to_char('5'_phred63), '5');
    EXPECT_EQ(to_char('6'_phred63), '6');
    EXPECT_EQ(to_char('7'_phred63), '7');
    EXPECT_EQ(to_char('8'_phred63), '8');
    EXPECT_EQ(to_char('9'_phred63), '9');
    EXPECT_EQ(to_char(':'_phred63), ':');
    EXPECT_EQ(to_char(';'_phred63), ';');
    EXPECT_EQ(to_char('<'_phred63), '<');
    EXPECT_EQ(to_char('='_phred63), '=');
    EXPECT_EQ(to_char('>'_phred63), '>');
    EXPECT_EQ(to_char('?'_phred63), '?');
    EXPECT_EQ(to_char('@'_phred63), '@');
    EXPECT_EQ(to_char('A'_phred63), 'A');
    EXPECT_EQ(to_char('B'_phred63), 'B');
    EXPECT_EQ(to_char('C'_phred63), 'C');
    EXPECT_EQ(to_char('D'_phred63), 'D');
    EXPECT_EQ(to_char('E'_phred63), 'E');
    EXPECT_EQ(to_char('F'_phred63), 'F');
    EXPECT_EQ(to_char('G'_phred63), 'G');
    EXPECT_EQ(to_char('H'_phred63), 'H');
    EXPECT_EQ(to_char('I'_phred63), 'I');
    EXPECT_EQ(to_char('J'_phred63), 'J');
    EXPECT_EQ(to_char('K'_phred63), 'K');
    EXPECT_EQ(to_char('L'_phred63), 'L');
    EXPECT_EQ(to_char('M'_phred63), 'M');
    EXPECT_EQ(to_char('N'_phred63), 'N');
    EXPECT_EQ(to_char('O'_phred63), 'O');
    EXPECT_EQ(to_char('P'_phred63), 'P');
    EXPECT_EQ(to_char('Q'_phred63), 'Q');
    EXPECT_EQ(to_char('R'_phred63), 'R');
    EXPECT_EQ(to_char('S'_phred63), 'S');
    EXPECT_EQ(to_char('T'_phred63), 'T');
    EXPECT_EQ(to_char('U'_phred63), 'U');
    EXPECT_EQ(to_char('V'_phred63), 'V');
    EXPECT_EQ(to_char('W'_phred63), 'W');
    EXPECT_EQ(to_char('X'_phred63), 'X');
    EXPECT_EQ(to_char('Y'_phred63), 'Y');
    EXPECT_EQ(to_char('Z'_phred63), 'Z');
    EXPECT_EQ(to_char('['_phred63), '[');
    EXPECT_EQ(to_char('\\'_phred63), '\\');
    EXPECT_EQ(to_char(']'_phred63), ']');
    EXPECT_EQ(to_char('^'_phred63), '^');
    EXPECT_EQ(to_char('_'_phred63), '_');
}

TEST(phred63, string_literal)
{
    std::vector<phred63> v;
    v.resize(5, '#'_phred63);
    EXPECT_EQ(v, "#####"_phred63);

    std::vector<phred63> w{'#'_phred63, '#'_phred63, '!'_phred63, '!'_phred63, '!'_phred63, '#'_phred63};
    EXPECT_EQ(w, "##!!!#"_phred63);
}

// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alphabet/quality/phred68legacy.hpp>

#include "../alphabet_constexpr_test_template.hpp"
#include "../alphabet_test_template.hpp"
#include "../semi_alphabet_constexpr_test_template.hpp"
#include "../semi_alphabet_test_template.hpp"
#include "phred_test_template.hpp"

using namespace seqan3;

INSTANTIATE_TYPED_TEST_SUITE_P(phred68legacy, alphabet_, phred68legacy, );
INSTANTIATE_TYPED_TEST_SUITE_P(phred68legacy, semi_alphabet_test, phred68legacy, );
INSTANTIATE_TYPED_TEST_SUITE_P(phred68legacy, alphabet_constexpr, phred68legacy, );
INSTANTIATE_TYPED_TEST_SUITE_P(phred68legacy, semi_alphabet_constexpr, phred68legacy, );
INSTANTIATE_TYPED_TEST_SUITE_P(phred68legacy, phred, phred68legacy, );

TEST(phred68legacy, char_literal)
{
    EXPECT_EQ(to_char(';'_phred68legacy), ';');
    EXPECT_EQ(to_char('<'_phred68legacy), '<');
    EXPECT_EQ(to_char('='_phred68legacy), '=');
    EXPECT_EQ(to_char('>'_phred68legacy), '>');
    EXPECT_EQ(to_char('?'_phred68legacy), '?');
    EXPECT_EQ(to_char('@'_phred68legacy), '@');
    EXPECT_EQ(to_char('A'_phred68legacy), 'A');
    EXPECT_EQ(to_char('B'_phred68legacy), 'B');
    EXPECT_EQ(to_char('C'_phred68legacy), 'C');
    EXPECT_EQ(to_char('D'_phred68legacy), 'D');
    EXPECT_EQ(to_char('E'_phred68legacy), 'E');
    EXPECT_EQ(to_char('F'_phred68legacy), 'F');
    EXPECT_EQ(to_char('G'_phred68legacy), 'G');
    EXPECT_EQ(to_char('H'_phred68legacy), 'H');
    EXPECT_EQ(to_char('I'_phred68legacy), 'I');
    EXPECT_EQ(to_char('J'_phred68legacy), 'J');
    EXPECT_EQ(to_char('K'_phred68legacy), 'K');
    EXPECT_EQ(to_char('L'_phred68legacy), 'L');
    EXPECT_EQ(to_char('M'_phred68legacy), 'M');
    EXPECT_EQ(to_char('N'_phred68legacy), 'N');
    EXPECT_EQ(to_char('O'_phred68legacy), 'O');
    EXPECT_EQ(to_char('P'_phred68legacy), 'P');
    EXPECT_EQ(to_char('Q'_phred68legacy), 'Q');
    EXPECT_EQ(to_char('R'_phred68legacy), 'R');
    EXPECT_EQ(to_char('S'_phred68legacy), 'S');
    EXPECT_EQ(to_char('T'_phred68legacy), 'T');
    EXPECT_EQ(to_char('U'_phred68legacy), 'U');
    EXPECT_EQ(to_char('V'_phred68legacy), 'V');
    EXPECT_EQ(to_char('W'_phred68legacy), 'W');
    EXPECT_EQ(to_char('X'_phred68legacy), 'X');
    EXPECT_EQ(to_char('Y'_phred68legacy), 'Y');
    EXPECT_EQ(to_char('Z'_phred68legacy), 'Z');
    EXPECT_EQ(to_char('['_phred68legacy), '[');
    EXPECT_EQ(to_char('\\'_phred68legacy), '\\');
    EXPECT_EQ(to_char(']'_phred68legacy), ']');
    EXPECT_EQ(to_char('^'_phred68legacy), '^');
    EXPECT_EQ(to_char('_'_phred68legacy), '_');
    EXPECT_EQ(to_char('`'_phred68legacy), '`');
    EXPECT_EQ(to_char('a'_phred68legacy), 'a');
    EXPECT_EQ(to_char('b'_phred68legacy), 'b');
    EXPECT_EQ(to_char('c'_phred68legacy), 'c');
    EXPECT_EQ(to_char('d'_phred68legacy), 'd');
    EXPECT_EQ(to_char('e'_phred68legacy), 'e');
    EXPECT_EQ(to_char('f'_phred68legacy), 'f');
    EXPECT_EQ(to_char('g'_phred68legacy), 'g');
    EXPECT_EQ(to_char('h'_phred68legacy), 'h');
    EXPECT_EQ(to_char('i'_phred68legacy), 'i');
    EXPECT_EQ(to_char('j'_phred68legacy), 'j');
    EXPECT_EQ(to_char('k'_phred68legacy), 'k');
    EXPECT_EQ(to_char('l'_phred68legacy), 'l');
    EXPECT_EQ(to_char('m'_phred68legacy), 'm');
    EXPECT_EQ(to_char('n'_phred68legacy), 'n');
    EXPECT_EQ(to_char('o'_phred68legacy), 'o');
    EXPECT_EQ(to_char('p'_phred68legacy), 'p');
    EXPECT_EQ(to_char('q'_phred68legacy), 'q');
    EXPECT_EQ(to_char('r'_phred68legacy), 'r');
    EXPECT_EQ(to_char('s'_phred68legacy), 's');
    EXPECT_EQ(to_char('t'_phred68legacy), 't');
    EXPECT_EQ(to_char('u'_phred68legacy), 'u');
    EXPECT_EQ(to_char('v'_phred68legacy), 'v');
    EXPECT_EQ(to_char('w'_phred68legacy), 'w');
    EXPECT_EQ(to_char('x'_phred68legacy), 'x');
    EXPECT_EQ(to_char('y'_phred68legacy), 'y');
    EXPECT_EQ(to_char('z'_phred68legacy), 'z');
    EXPECT_EQ(to_char('{'_phred68legacy), '{');
    EXPECT_EQ(to_char('|'_phred68legacy), '|');
    EXPECT_EQ(to_char('}'_phred68legacy), '}');
    EXPECT_EQ(to_char('~'_phred68legacy), '~');
}

TEST(phred68legacy, string_literal)
{
    std::vector<phred68legacy> v;
    v.resize(5, '#'_phred68legacy);
    EXPECT_EQ(v, "#####"_phred68legacy);

    std::vector<phred68legacy> w{'#'_phred68legacy, '#'_phred68legacy, '!'_phred68legacy, '!'_phred68legacy,
                                 '!'_phred68legacy, '#'_phred68legacy};
    EXPECT_EQ(w, "##!!!#"_phred68legacy);
}

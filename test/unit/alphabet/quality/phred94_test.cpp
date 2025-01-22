// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/alphabet/quality/phred94.hpp>
#include <seqan3/core/debug_stream/range.hpp>

#include "../alphabet_constexpr_test_template.hpp"
#include "../alphabet_test_template.hpp"
#include "../semi_alphabet_constexpr_test_template.hpp"
#include "../semi_alphabet_test_template.hpp"
#include "phred_test_template.hpp"

using seqan3::operator""_phred94;

INSTANTIATE_TYPED_TEST_SUITE_P(phred94, alphabet, seqan3::phred94, );
INSTANTIATE_TYPED_TEST_SUITE_P(phred94, semi_alphabet_test, seqan3::phred94, );
INSTANTIATE_TYPED_TEST_SUITE_P(phred94, alphabet_constexpr, seqan3::phred94, );
INSTANTIATE_TYPED_TEST_SUITE_P(phred94, semi_alphabet_constexpr, seqan3::phred94, );
INSTANTIATE_TYPED_TEST_SUITE_P(phred94, phred, seqan3::phred94, );

TEST(phred94, char_literal)
{
    EXPECT_EQ(seqan3::to_char('!'_phred94), '!');
    EXPECT_EQ(seqan3::to_char('"'_phred94), '"');
    EXPECT_EQ(seqan3::to_char('#'_phred94), '#');
    EXPECT_EQ(seqan3::to_char('$'_phred94), '$');
    EXPECT_EQ(seqan3::to_char('%'_phred94), '%');
    EXPECT_EQ(seqan3::to_char('&'_phred94), '&');
    EXPECT_EQ(seqan3::to_char('\''_phred94), '\'');
    EXPECT_EQ(seqan3::to_char('('_phred94), '(');
    EXPECT_EQ(seqan3::to_char(')'_phred94), ')');
    EXPECT_EQ(seqan3::to_char('*'_phred94), '*');
    EXPECT_EQ(seqan3::to_char('+'_phred94), '+');
    EXPECT_EQ(seqan3::to_char(','_phred94), ',');
    EXPECT_EQ(seqan3::to_char('-'_phred94), '-');
    EXPECT_EQ(seqan3::to_char('.'_phred94), '.');
    EXPECT_EQ(seqan3::to_char('/'_phred94), '/');
    EXPECT_EQ(seqan3::to_char('0'_phred94), '0');
    EXPECT_EQ(seqan3::to_char('1'_phred94), '1');
    EXPECT_EQ(seqan3::to_char('2'_phred94), '2');
    EXPECT_EQ(seqan3::to_char('3'_phred94), '3');
    EXPECT_EQ(seqan3::to_char('4'_phred94), '4');
    EXPECT_EQ(seqan3::to_char('5'_phred94), '5');
    EXPECT_EQ(seqan3::to_char('6'_phred94), '6');
    EXPECT_EQ(seqan3::to_char('7'_phred94), '7');
    EXPECT_EQ(seqan3::to_char('8'_phred94), '8');
    EXPECT_EQ(seqan3::to_char('9'_phred94), '9');
    EXPECT_EQ(seqan3::to_char(':'_phred94), ':');
    EXPECT_EQ(seqan3::to_char(';'_phred94), ';');
    EXPECT_EQ(seqan3::to_char('<'_phred94), '<');
    EXPECT_EQ(seqan3::to_char('='_phred94), '=');
    EXPECT_EQ(seqan3::to_char('>'_phred94), '>');
    EXPECT_EQ(seqan3::to_char('?'_phred94), '?');
    EXPECT_EQ(seqan3::to_char('@'_phred94), '@');
    EXPECT_EQ(seqan3::to_char('A'_phred94), 'A');
    EXPECT_EQ(seqan3::to_char('B'_phred94), 'B');
    EXPECT_EQ(seqan3::to_char('C'_phred94), 'C');
    EXPECT_EQ(seqan3::to_char('D'_phred94), 'D');
    EXPECT_EQ(seqan3::to_char('E'_phred94), 'E');
    EXPECT_EQ(seqan3::to_char('F'_phred94), 'F');
    EXPECT_EQ(seqan3::to_char('G'_phred94), 'G');
    EXPECT_EQ(seqan3::to_char('H'_phred94), 'H');
    EXPECT_EQ(seqan3::to_char('I'_phred94), 'I');
    EXPECT_EQ(seqan3::to_char('J'_phred94), 'J');
    EXPECT_EQ(seqan3::to_char('K'_phred94), 'K');
    EXPECT_EQ(seqan3::to_char('L'_phred94), 'L');
    EXPECT_EQ(seqan3::to_char('M'_phred94), 'M');
    EXPECT_EQ(seqan3::to_char('N'_phred94), 'N');
    EXPECT_EQ(seqan3::to_char('O'_phred94), 'O');
    EXPECT_EQ(seqan3::to_char('P'_phred94), 'P');
    EXPECT_EQ(seqan3::to_char('Q'_phred94), 'Q');
    EXPECT_EQ(seqan3::to_char('R'_phred94), 'R');
    EXPECT_EQ(seqan3::to_char('S'_phred94), 'S');
    EXPECT_EQ(seqan3::to_char('T'_phred94), 'T');
    EXPECT_EQ(seqan3::to_char('U'_phred94), 'U');
    EXPECT_EQ(seqan3::to_char('V'_phred94), 'V');
    EXPECT_EQ(seqan3::to_char('W'_phred94), 'W');
    EXPECT_EQ(seqan3::to_char('X'_phred94), 'X');
    EXPECT_EQ(seqan3::to_char('Y'_phred94), 'Y');
    EXPECT_EQ(seqan3::to_char('Z'_phred94), 'Z');
    EXPECT_EQ(seqan3::to_char('['_phred94), '[');
    EXPECT_EQ(seqan3::to_char('\\'_phred94), '\\');
    EXPECT_EQ(seqan3::to_char(']'_phred94), ']');
    EXPECT_EQ(seqan3::to_char('^'_phred94), '^');
    EXPECT_EQ(seqan3::to_char('_'_phred94), '_');
    EXPECT_EQ(seqan3::to_char('`'_phred94), '`');
    EXPECT_EQ(seqan3::to_char('a'_phred94), 'a');
    EXPECT_EQ(seqan3::to_char('b'_phred94), 'b');
    EXPECT_EQ(seqan3::to_char('c'_phred94), 'c');
    EXPECT_EQ(seqan3::to_char('d'_phred94), 'd');
    EXPECT_EQ(seqan3::to_char('e'_phred94), 'e');
    EXPECT_EQ(seqan3::to_char('f'_phred94), 'f');
    EXPECT_EQ(seqan3::to_char('g'_phred94), 'g');
    EXPECT_EQ(seqan3::to_char('h'_phred94), 'h');
    EXPECT_EQ(seqan3::to_char('i'_phred94), 'i');
    EXPECT_EQ(seqan3::to_char('j'_phred94), 'j');
    EXPECT_EQ(seqan3::to_char('k'_phred94), 'k');
    EXPECT_EQ(seqan3::to_char('l'_phred94), 'l');
    EXPECT_EQ(seqan3::to_char('m'_phred94), 'm');
    EXPECT_EQ(seqan3::to_char('n'_phred94), 'n');
    EXPECT_EQ(seqan3::to_char('o'_phred94), 'o');
    EXPECT_EQ(seqan3::to_char('p'_phred94), 'p');
    EXPECT_EQ(seqan3::to_char('q'_phred94), 'q');
    EXPECT_EQ(seqan3::to_char('r'_phred94), 'r');
    EXPECT_EQ(seqan3::to_char('s'_phred94), 's');
    EXPECT_EQ(seqan3::to_char('t'_phred94), 't');
    EXPECT_EQ(seqan3::to_char('u'_phred94), 'u');
    EXPECT_EQ(seqan3::to_char('v'_phred94), 'v');
    EXPECT_EQ(seqan3::to_char('w'_phred94), 'w');
    EXPECT_EQ(seqan3::to_char('x'_phred94), 'x');
    EXPECT_EQ(seqan3::to_char('y'_phred94), 'y');
    EXPECT_EQ(seqan3::to_char('z'_phred94), 'z');
    EXPECT_EQ(seqan3::to_char('{'_phred94), '{');
    EXPECT_EQ(seqan3::to_char('|'_phred94), '|');
    EXPECT_EQ(seqan3::to_char('}'_phred94), '}');
    EXPECT_EQ(seqan3::to_char('~'_phred94), '~');
}

TEST(phred94, string_literal)
{
    std::vector<seqan3::phred94> v;
    v.resize(5, '#'_phred94);
    EXPECT_EQ(v, "#####"_phred94);

    std::vector<seqan3::phred94> w{'#'_phred94, '#'_phred94, '!'_phred94, '!'_phred94, '!'_phred94, '#'_phred94};
    EXPECT_EQ(w, "##!!!#"_phred94);
}

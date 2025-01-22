// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/alphabet/quality/phred68solexa.hpp>
#include <seqan3/core/debug_stream/range.hpp>

#include "../alphabet_constexpr_test_template.hpp"
#include "../alphabet_test_template.hpp"
#include "../semi_alphabet_constexpr_test_template.hpp"
#include "../semi_alphabet_test_template.hpp"
#include "phred_test_template.hpp"

using seqan3::operator""_phred68solexa;

INSTANTIATE_TYPED_TEST_SUITE_P(phred68solexa, alphabet, seqan3::phred68solexa, );
INSTANTIATE_TYPED_TEST_SUITE_P(phred68solexa, semi_alphabet_test, seqan3::phred68solexa, );
INSTANTIATE_TYPED_TEST_SUITE_P(phred68solexa, alphabet_constexpr, seqan3::phred68solexa, );
INSTANTIATE_TYPED_TEST_SUITE_P(phred68solexa, semi_alphabet_constexpr, seqan3::phred68solexa, );
INSTANTIATE_TYPED_TEST_SUITE_P(phred68solexa, phred, seqan3::phred68solexa, );

TEST(phred68solexa, char_literal)
{
    EXPECT_EQ(seqan3::to_char(';'_phred68solexa), ';');
    EXPECT_EQ(seqan3::to_char('<'_phred68solexa), '<');
    EXPECT_EQ(seqan3::to_char('='_phred68solexa), '=');
    EXPECT_EQ(seqan3::to_char('>'_phred68solexa), '>');
    EXPECT_EQ(seqan3::to_char('?'_phred68solexa), '?');
    EXPECT_EQ(seqan3::to_char('@'_phred68solexa), '@');
    EXPECT_EQ(seqan3::to_char('A'_phred68solexa), 'A');
    EXPECT_EQ(seqan3::to_char('B'_phred68solexa), 'B');
    EXPECT_EQ(seqan3::to_char('C'_phred68solexa), 'C');
    EXPECT_EQ(seqan3::to_char('D'_phred68solexa), 'D');
    EXPECT_EQ(seqan3::to_char('E'_phred68solexa), 'E');
    EXPECT_EQ(seqan3::to_char('F'_phred68solexa), 'F');
    EXPECT_EQ(seqan3::to_char('G'_phred68solexa), 'G');
    EXPECT_EQ(seqan3::to_char('H'_phred68solexa), 'H');
    EXPECT_EQ(seqan3::to_char('I'_phred68solexa), 'I');
    EXPECT_EQ(seqan3::to_char('J'_phred68solexa), 'J');
    EXPECT_EQ(seqan3::to_char('K'_phred68solexa), 'K');
    EXPECT_EQ(seqan3::to_char('L'_phred68solexa), 'L');
    EXPECT_EQ(seqan3::to_char('M'_phred68solexa), 'M');
    EXPECT_EQ(seqan3::to_char('N'_phred68solexa), 'N');
    EXPECT_EQ(seqan3::to_char('O'_phred68solexa), 'O');
    EXPECT_EQ(seqan3::to_char('P'_phred68solexa), 'P');
    EXPECT_EQ(seqan3::to_char('Q'_phred68solexa), 'Q');
    EXPECT_EQ(seqan3::to_char('R'_phred68solexa), 'R');
    EXPECT_EQ(seqan3::to_char('S'_phred68solexa), 'S');
    EXPECT_EQ(seqan3::to_char('T'_phred68solexa), 'T');
    EXPECT_EQ(seqan3::to_char('U'_phred68solexa), 'U');
    EXPECT_EQ(seqan3::to_char('V'_phred68solexa), 'V');
    EXPECT_EQ(seqan3::to_char('W'_phred68solexa), 'W');
    EXPECT_EQ(seqan3::to_char('X'_phred68solexa), 'X');
    EXPECT_EQ(seqan3::to_char('Y'_phred68solexa), 'Y');
    EXPECT_EQ(seqan3::to_char('Z'_phred68solexa), 'Z');
    EXPECT_EQ(seqan3::to_char('['_phred68solexa), '[');
    EXPECT_EQ(seqan3::to_char('\\'_phred68solexa), '\\');
    EXPECT_EQ(seqan3::to_char(']'_phred68solexa), ']');
    EXPECT_EQ(seqan3::to_char('^'_phred68solexa), '^');
    EXPECT_EQ(seqan3::to_char('_'_phred68solexa), '_');
    EXPECT_EQ(seqan3::to_char('`'_phred68solexa), '`');
    EXPECT_EQ(seqan3::to_char('a'_phred68solexa), 'a');
    EXPECT_EQ(seqan3::to_char('b'_phred68solexa), 'b');
    EXPECT_EQ(seqan3::to_char('c'_phred68solexa), 'c');
    EXPECT_EQ(seqan3::to_char('d'_phred68solexa), 'd');
    EXPECT_EQ(seqan3::to_char('e'_phred68solexa), 'e');
    EXPECT_EQ(seqan3::to_char('f'_phred68solexa), 'f');
    EXPECT_EQ(seqan3::to_char('g'_phred68solexa), 'g');
    EXPECT_EQ(seqan3::to_char('h'_phred68solexa), 'h');
    EXPECT_EQ(seqan3::to_char('i'_phred68solexa), 'i');
    EXPECT_EQ(seqan3::to_char('j'_phred68solexa), 'j');
    EXPECT_EQ(seqan3::to_char('k'_phred68solexa), 'k');
    EXPECT_EQ(seqan3::to_char('l'_phred68solexa), 'l');
    EXPECT_EQ(seqan3::to_char('m'_phred68solexa), 'm');
    EXPECT_EQ(seqan3::to_char('n'_phred68solexa), 'n');
    EXPECT_EQ(seqan3::to_char('o'_phred68solexa), 'o');
    EXPECT_EQ(seqan3::to_char('p'_phred68solexa), 'p');
    EXPECT_EQ(seqan3::to_char('q'_phred68solexa), 'q');
    EXPECT_EQ(seqan3::to_char('r'_phred68solexa), 'r');
    EXPECT_EQ(seqan3::to_char('s'_phred68solexa), 's');
    EXPECT_EQ(seqan3::to_char('t'_phred68solexa), 't');
    EXPECT_EQ(seqan3::to_char('u'_phred68solexa), 'u');
    EXPECT_EQ(seqan3::to_char('v'_phred68solexa), 'v');
    EXPECT_EQ(seqan3::to_char('w'_phred68solexa), 'w');
    EXPECT_EQ(seqan3::to_char('x'_phred68solexa), 'x');
    EXPECT_EQ(seqan3::to_char('y'_phred68solexa), 'y');
    EXPECT_EQ(seqan3::to_char('z'_phred68solexa), 'z');
    EXPECT_EQ(seqan3::to_char('{'_phred68solexa), '{');
    EXPECT_EQ(seqan3::to_char('|'_phred68solexa), '|');
    EXPECT_EQ(seqan3::to_char('}'_phred68solexa), '}');
    EXPECT_EQ(seqan3::to_char('~'_phred68solexa), '~');
}

TEST(phred68solexa, string_literal)
{
    std::vector<seqan3::phred68solexa> v;
    v.resize(5, '#'_phred68solexa);
    EXPECT_EQ(v, "#####"_phred68solexa);

    std::vector<seqan3::phred68solexa> w{'#'_phred68solexa,
                                         '#'_phred68solexa,
                                         '!'_phred68solexa,
                                         '!'_phred68solexa,
                                         '!'_phred68solexa,
                                         '#'_phred68solexa};
    EXPECT_EQ(w, "##!!!#"_phred68solexa);
}

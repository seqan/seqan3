// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alphabet/quality/phred42.hpp>

#include "../alphabet_constexpr_test_template.hpp"
#include "../alphabet_test_template.hpp"
#include "../semi_alphabet_constexpr_test_template.hpp"
#include "../semi_alphabet_test_template.hpp"
#include "phred_test_template.hpp"

using namespace seqan3;

INSTANTIATE_TYPED_TEST_SUITE_P(phred42, alphabet_, phred42, );
INSTANTIATE_TYPED_TEST_SUITE_P(phred42, semi_alphabet_test, phred42, );
INSTANTIATE_TYPED_TEST_SUITE_P(phred42, alphabet_constexpr, phred42, );
INSTANTIATE_TYPED_TEST_SUITE_P(phred42, semi_alphabet_constexpr, phred42, );
INSTANTIATE_TYPED_TEST_SUITE_P(phred42, phred, phred42, );

TEST(phred42, char_literal)
{
    EXPECT_EQ(to_char('!'_phred42), '!');
    EXPECT_EQ(to_char('"'_phred42), '"');
    EXPECT_EQ(to_char('#'_phred42), '#');
    EXPECT_EQ(to_char('$'_phred42), '$');
    EXPECT_EQ(to_char('%'_phred42), '%');
    EXPECT_EQ(to_char('&'_phred42), '&');
    EXPECT_EQ(to_char('\''_phred42), '\'');
    EXPECT_EQ(to_char('('_phred42), '(');
    EXPECT_EQ(to_char(')'_phred42), ')');
    EXPECT_EQ(to_char('*'_phred42), '*');
    EXPECT_EQ(to_char('+'_phred42), '+');
    EXPECT_EQ(to_char(','_phred42), ',');
    EXPECT_EQ(to_char('-'_phred42), '-');
    EXPECT_EQ(to_char('.'_phred42), '.');
    EXPECT_EQ(to_char('/'_phred42), '/');
    EXPECT_EQ(to_char('0'_phred42), '0');
    EXPECT_EQ(to_char('1'_phred42), '1');
    EXPECT_EQ(to_char('2'_phred42), '2');
    EXPECT_EQ(to_char('3'_phred42), '3');
    EXPECT_EQ(to_char('4'_phred42), '4');
    EXPECT_EQ(to_char('5'_phred42), '5');
    EXPECT_EQ(to_char('6'_phred42), '6');
    EXPECT_EQ(to_char('7'_phred42), '7');
    EXPECT_EQ(to_char('8'_phred42), '8');
    EXPECT_EQ(to_char('9'_phred42), '9');
    EXPECT_EQ(to_char(':'_phred42), ':');
    EXPECT_EQ(to_char(';'_phred42), ';');
    EXPECT_EQ(to_char('<'_phred42), '<');
    EXPECT_EQ(to_char('='_phred42), '=');
    EXPECT_EQ(to_char('>'_phred42), '>');
    EXPECT_EQ(to_char('?'_phred42), '?');
    EXPECT_EQ(to_char('@'_phred42), '@');
    EXPECT_EQ(to_char('A'_phred42), 'A');
    EXPECT_EQ(to_char('B'_phred42), 'B');
    EXPECT_EQ(to_char('C'_phred42), 'C');
    EXPECT_EQ(to_char('D'_phred42), 'D');
    EXPECT_EQ(to_char('E'_phred42), 'E');
    EXPECT_EQ(to_char('F'_phred42), 'F');
    EXPECT_EQ(to_char('G'_phred42), 'G');
    EXPECT_EQ(to_char('H'_phred42), 'H');
    EXPECT_EQ(to_char('I'_phred42), 'I');
    EXPECT_EQ(to_char('J'_phred42), 'J');
}

TEST(phred42, string_literal)
{
    std::vector<phred42> v;
    v.resize(5, '#'_phred42);
    EXPECT_EQ(v, "#####"_phred42);

    std::vector<phred42> w{'#'_phred42, '#'_phred42, '!'_phred42, '!'_phred42, '!'_phred42, '#'_phred42};
    EXPECT_EQ(w, "##!!!#"_phred42);
}

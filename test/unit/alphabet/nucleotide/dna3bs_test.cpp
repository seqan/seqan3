// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <seqan3/alphabet/nucleotide/dna3bs.hpp>
#include <seqan3/core/char_operations/predicate.hpp>

#include "../alphabet_constexpr_test_template.hpp"
#include "../alphabet_test_template.hpp"
#include "../semi_alphabet_constexpr_test_template.hpp"
#include "../semi_alphabet_test_template.hpp"

INSTANTIATE_TYPED_TEST_SUITE_P(dna3bs, alphabet_, seqan3::dna3bs, );
INSTANTIATE_TYPED_TEST_SUITE_P(dna3bs, semi_alphabet_test, seqan3::dna3bs, );
INSTANTIATE_TYPED_TEST_SUITE_P(dna3bs, alphabet_constexpr, seqan3::dna3bs, );
INSTANTIATE_TYPED_TEST_SUITE_P(dna3bs, semi_alphabet_constexpr, seqan3::dna3bs, );

using seqan3::operator""_dna3bs;

TEST(dna3bs, concept_check)
{
    EXPECT_TRUE(seqan3::nucleotide_alphabet<seqan3::dna3bs>);
    EXPECT_TRUE(seqan3::nucleotide_alphabet<seqan3::dna3bs &>);
    EXPECT_TRUE(seqan3::nucleotide_alphabet<seqan3::dna3bs const>);
    EXPECT_TRUE(seqan3::nucleotide_alphabet<seqan3::dna3bs const &>);
}

TEST(dna3bs, global_complement)
{
    EXPECT_EQ(seqan3::complement(seqan3::dna3bs{}.assign_char('A')), seqan3::dna3bs{}.assign_char('T'));
    EXPECT_EQ(seqan3::complement(seqan3::dna3bs{}.assign_char('C')), seqan3::dna3bs{}.assign_char('A'));
    EXPECT_EQ(seqan3::complement(seqan3::dna3bs{}.assign_char('G')), seqan3::dna3bs{}.assign_char('T'));
    EXPECT_EQ(seqan3::complement(seqan3::dna3bs{}.assign_char('T')), seqan3::dna3bs{}.assign_char('A'));
}

TEST(dna3bs, to_char_assign_char)
{
    EXPECT_EQ(seqan3::to_char(seqan3::dna3bs{}.assign_char('A')), 'A');
    EXPECT_EQ(seqan3::to_char(seqan3::dna3bs{}.assign_char('C')), 'T');
    EXPECT_EQ(seqan3::to_char(seqan3::dna3bs{}.assign_char('G')), 'G');

    EXPECT_EQ(seqan3::to_char(seqan3::dna3bs{}.assign_char('U')), 'T');
    EXPECT_EQ(seqan3::to_char(seqan3::dna3bs{}.assign_char('T')), 'T');

    EXPECT_EQ(seqan3::to_char(seqan3::dna3bs{}.assign_char('R')), 'A');
    EXPECT_EQ(seqan3::to_char(seqan3::dna3bs{}.assign_char('Y')), 'T');
    EXPECT_EQ(seqan3::to_char(seqan3::dna3bs{}.assign_char('S')), 'T');
    EXPECT_EQ(seqan3::to_char(seqan3::dna3bs{}.assign_char('W')), 'A');
    EXPECT_EQ(seqan3::to_char(seqan3::dna3bs{}.assign_char('K')), 'G');
    EXPECT_EQ(seqan3::to_char(seqan3::dna3bs{}.assign_char('M')), 'A');
    EXPECT_EQ(seqan3::to_char(seqan3::dna3bs{}.assign_char('B')), 'T');
    EXPECT_EQ(seqan3::to_char(seqan3::dna3bs{}.assign_char('D')), 'A');
    EXPECT_EQ(seqan3::to_char(seqan3::dna3bs{}.assign_char('H')), 'A');
    EXPECT_EQ(seqan3::to_char(seqan3::dna3bs{}.assign_char('V')), 'A');

    EXPECT_EQ(seqan3::to_char(seqan3::dna3bs{}.assign_char('N')), 'A');
    EXPECT_EQ(seqan3::to_char(seqan3::dna3bs{}.assign_char('!')), 'A');
}

TEST(dna3bs, char_literal)
{
    EXPECT_EQ(seqan3::to_char('A'_dna3bs), 'A');
    EXPECT_EQ(seqan3::to_char('C'_dna3bs), 'T');
    EXPECT_EQ(seqan3::to_char('G'_dna3bs), 'G');

    EXPECT_EQ(seqan3::to_char('U'_dna3bs), 'T');
    EXPECT_EQ(seqan3::to_char('T'_dna3bs), 'T');

    EXPECT_EQ(seqan3::to_char('R'_dna3bs), 'A');
    EXPECT_EQ(seqan3::to_char('Y'_dna3bs), 'T');
    EXPECT_EQ(seqan3::to_char('S'_dna3bs), 'T');
    EXPECT_EQ(seqan3::to_char('W'_dna3bs), 'A');
    EXPECT_EQ(seqan3::to_char('K'_dna3bs), 'G');
    EXPECT_EQ(seqan3::to_char('M'_dna3bs), 'A');
    EXPECT_EQ(seqan3::to_char('B'_dna3bs), 'T');
    EXPECT_EQ(seqan3::to_char('D'_dna3bs), 'A');
    EXPECT_EQ(seqan3::to_char('H'_dna3bs), 'A');
    EXPECT_EQ(seqan3::to_char('V'_dna3bs), 'A');

    EXPECT_EQ(seqan3::to_char('N'_dna3bs), 'A');
    EXPECT_EQ(seqan3::to_char('!'_dna3bs), 'A');
}

TEST(dna3bs, string_literal)
{
    seqan3::dna3bs_vector v;
    v.resize(5, 'A'_dna3bs);
    EXPECT_EQ(v, "AAAAA"_dna3bs);

    std::vector<seqan3::dna3bs> w{'A'_dna3bs, 'C'_dna3bs, 'G'_dna3bs, 'T'_dna3bs, 'U'_dna3bs, 'N'_dna3bs};
    EXPECT_EQ(w, "ATGTTA"_dna3bs);
}

TEST(dna3bs, char_is_valid)
{
    constexpr auto validator = seqan3::is_char<'A'> || seqan3::is_char<'G'> || seqan3::is_char<'T'>
                            || seqan3::is_char<'U'> || seqan3::is_char<'a'> || seqan3::is_char<'g'>
                            || seqan3::is_char<'t'> || seqan3::is_char<'u'>;
    for (char c : std::views::iota(std::numeric_limits<char>::min(), std::numeric_limits<char>::max()))
        EXPECT_EQ(seqan3::dna3bs::char_is_valid(c), validator(c));
}

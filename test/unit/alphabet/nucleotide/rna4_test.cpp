// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <seqan3/alphabet/nucleotide/rna4.hpp>
#include <seqan3/core/debug_stream/range.hpp>
#include <seqan3/utility/char_operations/predicate.hpp>

#include "../alphabet_constexpr_test_template.hpp"
#include "../alphabet_test_template.hpp"
#include "../semi_alphabet_constexpr_test_template.hpp"
#include "../semi_alphabet_test_template.hpp"
#include "nucleotide_test_template.hpp"

using seqan3::operator""_rna4;

INSTANTIATE_TYPED_TEST_SUITE_P(rna4, alphabet, seqan3::rna4, );
INSTANTIATE_TYPED_TEST_SUITE_P(rna4, semi_alphabet_test, seqan3::rna4, );
INSTANTIATE_TYPED_TEST_SUITE_P(rna4, alphabet_constexpr, seqan3::rna4, );
INSTANTIATE_TYPED_TEST_SUITE_P(rna4, semi_alphabet_constexpr, seqan3::rna4, );
INSTANTIATE_TYPED_TEST_SUITE_P(rna4, nucleotide, seqan3::rna4, );

TEST(rna4, to_char_assign_char)
{
    EXPECT_EQ(seqan3::to_char(seqan3::rna4{}.assign_char('A')), 'A');
    EXPECT_EQ(seqan3::to_char(seqan3::rna4{}.assign_char('C')), 'C');
    EXPECT_EQ(seqan3::to_char(seqan3::rna4{}.assign_char('G')), 'G');

    EXPECT_EQ(seqan3::to_char(seqan3::rna4{}.assign_char('U')), 'U');
    EXPECT_EQ(seqan3::to_char(seqan3::rna4{}.assign_char('T')), 'U');

    EXPECT_EQ(seqan3::to_char(seqan3::rna4{}.assign_char('R')), 'A');
    EXPECT_EQ(seqan3::to_char(seqan3::rna4{}.assign_char('Y')), 'C');
    EXPECT_EQ(seqan3::to_char(seqan3::rna4{}.assign_char('S')), 'C');
    EXPECT_EQ(seqan3::to_char(seqan3::rna4{}.assign_char('W')), 'A');
    EXPECT_EQ(seqan3::to_char(seqan3::rna4{}.assign_char('K')), 'G');
    EXPECT_EQ(seqan3::to_char(seqan3::rna4{}.assign_char('M')), 'A');
    EXPECT_EQ(seqan3::to_char(seqan3::rna4{}.assign_char('B')), 'C');
    EXPECT_EQ(seqan3::to_char(seqan3::rna4{}.assign_char('D')), 'A');
    EXPECT_EQ(seqan3::to_char(seqan3::rna4{}.assign_char('H')), 'A');
    EXPECT_EQ(seqan3::to_char(seqan3::rna4{}.assign_char('V')), 'A');

    EXPECT_EQ(seqan3::to_char(seqan3::rna4{}.assign_char('N')), 'A');
    EXPECT_EQ(seqan3::to_char(seqan3::rna4{}.assign_char('!')), 'A');
}

TEST(rna4, char_literal)
{
    EXPECT_EQ(seqan3::to_char('A'_rna4), 'A');
    EXPECT_EQ(seqan3::to_char('C'_rna4), 'C');
    EXPECT_EQ(seqan3::to_char('G'_rna4), 'G');

    EXPECT_EQ(seqan3::to_char('U'_rna4), 'U');
    EXPECT_EQ(seqan3::to_char('T'_rna4), 'U');

    EXPECT_EQ(seqan3::to_char('R'_rna4), 'A');
    EXPECT_EQ(seqan3::to_char('Y'_rna4), 'C');
    EXPECT_EQ(seqan3::to_char('S'_rna4), 'C');
    EXPECT_EQ(seqan3::to_char('W'_rna4), 'A');
    EXPECT_EQ(seqan3::to_char('K'_rna4), 'G');
    EXPECT_EQ(seqan3::to_char('M'_rna4), 'A');
    EXPECT_EQ(seqan3::to_char('B'_rna4), 'C');
    EXPECT_EQ(seqan3::to_char('D'_rna4), 'A');
    EXPECT_EQ(seqan3::to_char('H'_rna4), 'A');
    EXPECT_EQ(seqan3::to_char('V'_rna4), 'A');

    EXPECT_EQ(seqan3::to_char('N'_rna4), 'A');
    EXPECT_EQ(seqan3::to_char('!'_rna4), 'A');
}

TEST(rna4, string_literal)
{
    seqan3::rna4_vector v;
    v.resize(5, 'A'_rna4);
    EXPECT_EQ(v, "AAAAA"_rna4);

    std::vector<seqan3::rna4> w{'A'_rna4, 'C'_rna4, 'G'_rna4, 'T'_rna4, 'U'_rna4, 'N'_rna4};
    EXPECT_EQ(w, "ACGUUA"_rna4);
}

TEST(rna4, char_is_valid)
{
    constexpr auto validator = seqan3::is_char<'A'> || seqan3::is_char<'C'> || seqan3::is_char<'G'>
                            || seqan3::is_char<'T'> || seqan3::is_char<'U'> || seqan3::is_char<'a'>
                            || seqan3::is_char<'c'> || seqan3::is_char<'g'> || seqan3::is_char<'t'>
                            || seqan3::is_char<'u'>;
    for (char c : std::views::iota(std::numeric_limits<char>::min(), std::numeric_limits<char>::max()))
        EXPECT_EQ(seqan3::rna4::char_is_valid(c), validator(c));
}

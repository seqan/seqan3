// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <seqan3/alphabet/nucleotide/rna5.hpp>
#include <seqan3/core/debug_stream/range.hpp>
#include <seqan3/utility/char_operations/predicate.hpp>

#include "../alphabet_constexpr_test_template.hpp"
#include "../alphabet_test_template.hpp"
#include "../semi_alphabet_constexpr_test_template.hpp"
#include "../semi_alphabet_test_template.hpp"
#include "nucleotide_test_template.hpp"

using seqan3::operator""_rna5;

INSTANTIATE_TYPED_TEST_SUITE_P(rna5, alphabet, seqan3::rna5, );
INSTANTIATE_TYPED_TEST_SUITE_P(rna5, semi_alphabet_test, seqan3::rna5, );
INSTANTIATE_TYPED_TEST_SUITE_P(rna5, alphabet_constexpr, seqan3::rna5, );
INSTANTIATE_TYPED_TEST_SUITE_P(rna5, semi_alphabet_constexpr, seqan3::rna5, );
INSTANTIATE_TYPED_TEST_SUITE_P(rna5, nucleotide, seqan3::rna5, );

TEST(rna5, to_char_assign_char)
{
    EXPECT_EQ(seqan3::to_char(seqan3::rna5{}.assign_char('A')), 'A');
    EXPECT_EQ(seqan3::to_char(seqan3::rna5{}.assign_char('C')), 'C');
    EXPECT_EQ(seqan3::to_char(seqan3::rna5{}.assign_char('G')), 'G');

    EXPECT_EQ(seqan3::to_char(seqan3::rna5{}.assign_char('U')), 'U');
    EXPECT_EQ(seqan3::to_char(seqan3::rna5{}.assign_char('T')), 'U');

    EXPECT_EQ(seqan3::to_char(seqan3::rna5{}.assign_char('R')), 'N');
    EXPECT_EQ(seqan3::to_char(seqan3::rna5{}.assign_char('Y')), 'N');
    EXPECT_EQ(seqan3::to_char(seqan3::rna5{}.assign_char('S')), 'N');
    EXPECT_EQ(seqan3::to_char(seqan3::rna5{}.assign_char('W')), 'N');
    EXPECT_EQ(seqan3::to_char(seqan3::rna5{}.assign_char('K')), 'N');
    EXPECT_EQ(seqan3::to_char(seqan3::rna5{}.assign_char('M')), 'N');
    EXPECT_EQ(seqan3::to_char(seqan3::rna5{}.assign_char('B')), 'N');
    EXPECT_EQ(seqan3::to_char(seqan3::rna5{}.assign_char('D')), 'N');
    EXPECT_EQ(seqan3::to_char(seqan3::rna5{}.assign_char('H')), 'N');
    EXPECT_EQ(seqan3::to_char(seqan3::rna5{}.assign_char('V')), 'N');

    EXPECT_EQ(seqan3::to_char(seqan3::rna5{}.assign_char('N')), 'N');
    EXPECT_EQ(seqan3::to_char(seqan3::rna5{}.assign_char('!')), 'N');
}

TEST(rna5, char_literal)
{
    EXPECT_EQ(seqan3::to_char('A'_rna5), 'A');
    EXPECT_EQ(seqan3::to_char('C'_rna5), 'C');
    EXPECT_EQ(seqan3::to_char('G'_rna5), 'G');

    EXPECT_EQ(seqan3::to_char('U'_rna5), 'U');
    EXPECT_EQ(seqan3::to_char('T'_rna5), 'U');

    EXPECT_EQ(seqan3::to_char('R'_rna5), 'N');
    EXPECT_EQ(seqan3::to_char('Y'_rna5), 'N');
    EXPECT_EQ(seqan3::to_char('S'_rna5), 'N');
    EXPECT_EQ(seqan3::to_char('W'_rna5), 'N');
    EXPECT_EQ(seqan3::to_char('K'_rna5), 'N');
    EXPECT_EQ(seqan3::to_char('M'_rna5), 'N');
    EXPECT_EQ(seqan3::to_char('B'_rna5), 'N');
    EXPECT_EQ(seqan3::to_char('D'_rna5), 'N');
    EXPECT_EQ(seqan3::to_char('H'_rna5), 'N');
    EXPECT_EQ(seqan3::to_char('V'_rna5), 'N');

    EXPECT_EQ(seqan3::to_char('N'_rna5), 'N');
    EXPECT_EQ(seqan3::to_char('!'_rna5), 'N');
}

TEST(rna5, string_literal)
{
    seqan3::rna5_vector v;
    v.resize(5, 'A'_rna5);
    EXPECT_EQ(v, "AAAAA"_rna5);

    std::vector<seqan3::rna5> w{'A'_rna5, 'C'_rna5, 'G'_rna5, 'T'_rna5, 'U'_rna5, 'N'_rna5};
    EXPECT_EQ(w, "ACGUUN"_rna5);
}

TEST(rna5, char_is_valid)
{
    constexpr auto validator = seqan3::is_char<'A'> || seqan3::is_char<'C'> || seqan3::is_char<'G'>
                            || seqan3::is_char<'T'> || seqan3::is_char<'U'> || seqan3::is_char<'a'>
                            || seqan3::is_char<'c'> || seqan3::is_char<'g'> || seqan3::is_char<'t'>
                            || seqan3::is_char<'u'> || seqan3::is_char<'N'> || seqan3::is_char<'n'>;
    for (char c : std::views::iota(std::numeric_limits<char>::min(), std::numeric_limits<char>::max()))
        EXPECT_EQ(seqan3::rna5::char_is_valid(c), validator(c));
}

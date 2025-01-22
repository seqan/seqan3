// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <seqan3/alphabet/nucleotide/dna16sam.hpp>
#include <seqan3/core/debug_stream/range.hpp>
#include <seqan3/utility/char_operations/predicate.hpp>

#include "../alphabet_constexpr_test_template.hpp"
#include "../alphabet_test_template.hpp"
#include "../semi_alphabet_constexpr_test_template.hpp"
#include "../semi_alphabet_test_template.hpp"

using seqan3::operator""_dna16sam;

// ------------------------------------------------------------------
// dna16sam alphabet
// ------------------------------------------------------------------

INSTANTIATE_TYPED_TEST_SUITE_P(dna16sam, alphabet, seqan3::dna16sam, );
INSTANTIATE_TYPED_TEST_SUITE_P(dna16sam, semi_alphabet_test, seqan3::dna16sam, );
INSTANTIATE_TYPED_TEST_SUITE_P(dna16sam, alphabet_constexpr, seqan3::dna16sam, );
INSTANTIATE_TYPED_TEST_SUITE_P(dna16sam, semi_alphabet_constexpr, seqan3::dna16sam, );

// nucleotide test: (because the complement is not bijective for dna16sam we need to test it manually)
TEST(dna16sam, nucleotide)
{
    EXPECT_TRUE(seqan3::nucleotide_alphabet<seqan3::dna16sam>);
    EXPECT_TRUE(seqan3::nucleotide_alphabet<seqan3::dna16sam &>);

    EXPECT_EQ(seqan3::complement('='_dna16sam), 'N'_dna16sam);
    EXPECT_EQ(seqan3::complement('A'_dna16sam), 'T'_dna16sam);
    EXPECT_EQ(seqan3::complement('C'_dna16sam), 'G'_dna16sam);
    EXPECT_EQ(seqan3::complement('M'_dna16sam), 'K'_dna16sam);
    EXPECT_EQ(seqan3::complement('G'_dna16sam), 'C'_dna16sam);
    EXPECT_EQ(seqan3::complement('R'_dna16sam), 'Y'_dna16sam);
    EXPECT_EQ(seqan3::complement('S'_dna16sam), 'S'_dna16sam);
    EXPECT_EQ(seqan3::complement('V'_dna16sam), 'B'_dna16sam);
    EXPECT_EQ(seqan3::complement('T'_dna16sam), 'A'_dna16sam);
    EXPECT_EQ(seqan3::complement('W'_dna16sam), 'W'_dna16sam);
    EXPECT_EQ(seqan3::complement('Y'_dna16sam), 'R'_dna16sam);
    EXPECT_EQ(seqan3::complement('H'_dna16sam), 'D'_dna16sam);
    EXPECT_EQ(seqan3::complement('K'_dna16sam), 'M'_dna16sam);
    EXPECT_EQ(seqan3::complement('D'_dna16sam), 'H'_dna16sam);
    EXPECT_EQ(seqan3::complement('B'_dna16sam), 'V'_dna16sam);
    EXPECT_EQ(seqan3::complement('N'_dna16sam), 'N'_dna16sam);
}

TEST(dna16sam, to_char_assign_char)
{
    using rank_t = seqan3::alphabet_rank_t<seqan3::dna16sam>;
    for (rank_t rank = 0; rank < seqan3::alphabet_size<seqan3::dna16sam>; ++rank)
    {
        char chr = seqan3::to_char(seqan3::assign_rank_to(rank, seqan3::dna16sam{}));
        EXPECT_EQ(seqan3::to_char(seqan3::dna16sam{}.assign_char(chr)), chr);
    }

    EXPECT_EQ(seqan3::to_char(seqan3::dna16sam{}.assign_char('a')), 'A');
    EXPECT_EQ(seqan3::to_char(seqan3::dna16sam{}.assign_char('c')), 'C');
    EXPECT_EQ(seqan3::to_char(seqan3::dna16sam{}.assign_char('g')), 'G');
    EXPECT_EQ(seqan3::to_char(seqan3::dna16sam{}.assign_char('t')), 'T');

    EXPECT_EQ(seqan3::to_char(seqan3::dna16sam{}.assign_char('U')), 'T');
    EXPECT_EQ(seqan3::to_char(seqan3::dna16sam{}.assign_char('!')), 'N');
}

TEST(dna16sam, char_literal)
{
    EXPECT_EQ(seqan3::to_char('A'_dna16sam), 'A');
    EXPECT_EQ(seqan3::to_char('C'_dna16sam), 'C');
    EXPECT_EQ(seqan3::to_char('G'_dna16sam), 'G');

    EXPECT_EQ(seqan3::to_char('U'_dna16sam), 'T');
    EXPECT_EQ(seqan3::to_char('T'_dna16sam), 'T');

    EXPECT_EQ(seqan3::to_char('R'_dna16sam), 'R');
    EXPECT_EQ(seqan3::to_char('Y'_dna16sam), 'Y');
    EXPECT_EQ(seqan3::to_char('S'_dna16sam), 'S');
    EXPECT_EQ(seqan3::to_char('W'_dna16sam), 'W');
    EXPECT_EQ(seqan3::to_char('K'_dna16sam), 'K');
    EXPECT_EQ(seqan3::to_char('M'_dna16sam), 'M');
    EXPECT_EQ(seqan3::to_char('B'_dna16sam), 'B');
    EXPECT_EQ(seqan3::to_char('D'_dna16sam), 'D');
    EXPECT_EQ(seqan3::to_char('H'_dna16sam), 'H');
    EXPECT_EQ(seqan3::to_char('V'_dna16sam), 'V');

    EXPECT_EQ(seqan3::to_char('='_dna16sam), '=');

    EXPECT_EQ(seqan3::to_char('N'_dna16sam), 'N');
    EXPECT_EQ(seqan3::to_char('!'_dna16sam), 'N');
}

TEST(dna16sam, string_literal)
{
    seqan3::dna16sam_vector v;
    v.resize(5, 'A'_dna16sam);
    EXPECT_EQ(v, "AAAAA"_dna16sam);

    std::vector<seqan3::dna16sam> w{'A'_dna16sam, '='_dna16sam, 'G'_dna16sam, 'T'_dna16sam, 'U'_dna16sam, 'N'_dna16sam};
    EXPECT_EQ(w, "A=GTTN"_dna16sam);
}

TEST(dna16sam, char_is_valid)
{
    constexpr auto validator =
        seqan3::is_char<'A'> || seqan3::is_char<'C'> || seqan3::is_char<'G'> || seqan3::is_char<'T'>
        || seqan3::is_char<'U'> || seqan3::is_char<'a'> || seqan3::is_char<'c'> || seqan3::is_char<'g'>
        || seqan3::is_char<'t'> || seqan3::is_char<'u'> || seqan3::is_char<'N'> || seqan3::is_char<'n'>
        || seqan3::is_char<'R'> || seqan3::is_char<'Y'> || seqan3::is_char<'S'> || seqan3::is_char<'W'>
        || seqan3::is_char<'K'> || seqan3::is_char<'M'> || seqan3::is_char<'B'> || seqan3::is_char<'D'>
        || seqan3::is_char<'H'> || seqan3::is_char<'V'> || seqan3::is_char<'r'> || seqan3::is_char<'y'>
        || seqan3::is_char<'s'> || seqan3::is_char<'w'> || seqan3::is_char<'k'> || seqan3::is_char<'m'>
        || seqan3::is_char<'b'> || seqan3::is_char<'d'> || seqan3::is_char<'h'> || seqan3::is_char<'v'>
        || seqan3::is_char<'='>;

    for (char c : std::views::iota(std::numeric_limits<char>::min(), std::numeric_limits<char>::max()))
        EXPECT_EQ(seqan3::dna16sam::char_is_valid(c), validator(c));
}

// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <seqan3/alphabet/structure/concept.hpp>
#include <seqan3/alphabet/structure/wuss.hpp>
#include <seqan3/range/views/zip.hpp>

#include "../alphabet_constexpr_test_template.hpp"
#include "../alphabet_test_template.hpp"
#include "../semi_alphabet_constexpr_test_template.hpp"
#include "../semi_alphabet_test_template.hpp"

using seqan3::operator""_wuss51;

INSTANTIATE_TYPED_TEST_SUITE_P(wuss51, alphabet_, seqan3::wuss51, );
INSTANTIATE_TYPED_TEST_SUITE_P(wuss51, semi_alphabet_test, seqan3::wuss51, );
INSTANTIATE_TYPED_TEST_SUITE_P(wuss51, alphabet_constexpr, seqan3::wuss51, );
INSTANTIATE_TYPED_TEST_SUITE_P(wuss51, semi_alphabet_constexpr, seqan3::wuss51, );
INSTANTIATE_TYPED_TEST_SUITE_P(wuss15, alphabet_, seqan3::wuss<15>, );
INSTANTIATE_TYPED_TEST_SUITE_P(wuss15, semi_alphabet_test, seqan3::wuss<15>, );
INSTANTIATE_TYPED_TEST_SUITE_P(wuss15, alphabet_constexpr, seqan3::wuss<15>, );
INSTANTIATE_TYPED_TEST_SUITE_P(wuss15, semi_alphabet_constexpr, seqan3::wuss<15>, );
INSTANTIATE_TYPED_TEST_SUITE_P(wuss67, alphabet_, seqan3::wuss<67>, );
INSTANTIATE_TYPED_TEST_SUITE_P(wuss67, semi_alphabet_test, seqan3::wuss<67>, );
INSTANTIATE_TYPED_TEST_SUITE_P(wuss67, alphabet_constexpr, seqan3::wuss<67>, );
INSTANTIATE_TYPED_TEST_SUITE_P(wuss67, semi_alphabet_constexpr, seqan3::wuss<67>, );

// assign_char functions
TEST(wuss, assign_char)
{
    std::vector<char> input
    {
        '.', '(', ')',
        ':', ',', '-', '_', '~', ';',
        '<', '>', '[', ']', '{', '}',
        'H', 'B', 'E', 'G', 'I', 'T', 'S'
    };

    std::vector<seqan3::wuss51> cmp
    {
        '.'_wuss51, '('_wuss51, ')'_wuss51,
        ':'_wuss51, ','_wuss51, '-'_wuss51, '_'_wuss51, '~'_wuss51, ';'_wuss51,
        '<'_wuss51, '>'_wuss51, '['_wuss51, ']'_wuss51, '{'_wuss51, '}'_wuss51,
        'H'_wuss51, 'B'_wuss51, 'E'_wuss51, 'G'_wuss51, 'I'_wuss51, 'T'_wuss51, 'S'_wuss51
    };

    for (auto [ ch, cm ] : seqan3::views::zip(input, cmp))
        EXPECT_EQ((seqan3::assign_char_to(ch, seqan3::wuss51{})), cm);
}

// to_char functions
TEST(wuss, to_char)
{
    EXPECT_EQ(seqan3::to_char('.'_wuss51), '.');
    EXPECT_EQ(seqan3::to_char(':'_wuss51), ':');
    EXPECT_EQ(seqan3::to_char(','_wuss51), ',');
    EXPECT_EQ(seqan3::to_char('-'_wuss51), '-');
    EXPECT_EQ(seqan3::to_char('_'_wuss51), '_');
    EXPECT_EQ(seqan3::to_char('~'_wuss51), '~');
    EXPECT_EQ(seqan3::to_char(';'_wuss51), ';');
    EXPECT_EQ(seqan3::to_char('<'_wuss51), '<');
    EXPECT_EQ(seqan3::to_char('>'_wuss51), '>');
    EXPECT_EQ(seqan3::to_char('('_wuss51), '(');
    EXPECT_EQ(seqan3::to_char(')'_wuss51), ')');
    EXPECT_EQ(seqan3::to_char('['_wuss51), '[');
    EXPECT_EQ(seqan3::to_char(']'_wuss51), ']');
    EXPECT_EQ(seqan3::to_char('{'_wuss51), '{');
    EXPECT_EQ(seqan3::to_char('}'_wuss51), '}');
}

// concepts
TEST(wuss, concept_check)
{
    EXPECT_TRUE(seqan3::rna_structure_alphabet<seqan3::wuss51>);
    EXPECT_TRUE(seqan3::rna_structure_alphabet<seqan3::wuss51 &>);
    EXPECT_TRUE(seqan3::rna_structure_alphabet<seqan3::wuss51 const>);
    EXPECT_TRUE(seqan3::rna_structure_alphabet<seqan3::wuss51 const &>);
    EXPECT_NE(seqan3::max_pseudoknot_depth<seqan3::wuss51>, 0);

    EXPECT_TRUE(seqan3::rna_structure_alphabet<seqan3::wuss<>>);  // same as wuss51
    EXPECT_TRUE(seqan3::rna_structure_alphabet<seqan3::wuss<> &>);
    EXPECT_TRUE(seqan3::rna_structure_alphabet<seqan3::wuss<> const>);
    EXPECT_TRUE(seqan3::rna_structure_alphabet<seqan3::wuss<> const &>);
    EXPECT_NE(seqan3::max_pseudoknot_depth<seqan3::wuss<>>, 0);

    EXPECT_TRUE(seqan3::rna_structure_alphabet<seqan3::wuss<67>>);
    EXPECT_TRUE(seqan3::rna_structure_alphabet<seqan3::wuss<67> &>);
    EXPECT_TRUE(seqan3::rna_structure_alphabet<seqan3::wuss<67> const>);
    EXPECT_TRUE(seqan3::rna_structure_alphabet<seqan3::wuss<67> const &>);
    EXPECT_NE(seqan3::max_pseudoknot_depth<seqan3::wuss<67>>, 0);

}

TEST(wuss, literals)
{
    std::vector<seqan3::wuss51> vec1;
    vec1.resize(5, '<'_wuss51);
    EXPECT_EQ(vec1, "<<<<<"_wuss51);

    std::vector<seqan3::wuss51> vec2{'.'_wuss51, '<'_wuss51, '<'_wuss51, '>'_wuss51, '>'_wuss51, '.'_wuss51};
    EXPECT_EQ(vec2, ".<<>>."_wuss51);
}

TEST(wuss, rna_structure_properties)
{
    EXPECT_EQ(seqan3::wuss51::max_pseudoknot_depth, 22);
    std::vector<seqan3::wuss51> vec = ".:,-_~;<>()[]{}AaBbCcDd"_wuss51;
    for (unsigned idx = 0; idx <= 6; ++idx)
    {
        EXPECT_TRUE(vec[idx].is_unpaired());
        EXPECT_FALSE(vec[idx].is_pair_open());
        EXPECT_FALSE(vec[idx].is_pair_close());
        EXPECT_FALSE(vec[idx].pseudoknot_id());
    }
    for (unsigned idx = 7; idx <= 21; idx+=2)
    {
        EXPECT_TRUE(vec[idx].is_pair_open());
        EXPECT_FALSE(vec[idx].is_unpaired());
        EXPECT_FALSE(vec[idx].is_pair_close());
        EXPECT_TRUE(vec[idx].pseudoknot_id());
        EXPECT_EQ(vec[idx].pseudoknot_id(), std::optional<uint8_t>{(idx - 7) / 2});
    }
    for (unsigned idx = 8; idx <= 22; idx+=2)
    {
        EXPECT_TRUE(vec[idx].is_pair_close());
        EXPECT_FALSE(vec[idx].is_unpaired());
        EXPECT_FALSE(vec[idx].is_pair_open());
        EXPECT_TRUE(vec[idx].pseudoknot_id());
        EXPECT_EQ(vec[idx].pseudoknot_id(), std::optional<uint8_t>{(idx - 8) / 2});
    }
}

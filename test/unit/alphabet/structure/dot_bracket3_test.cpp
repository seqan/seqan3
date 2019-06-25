// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <range/v3/view/zip.hpp>

#include <seqan3/alphabet/structure/concept.hpp>
#include <seqan3/alphabet/structure/dot_bracket3.hpp>

#include "../alphabet_test_template.hpp"
#include "../alphabet_constexpr_test_template.hpp"

INSTANTIATE_TYPED_TEST_CASE_P(dot_bracket3, alphabet, dot_bracket3);
INSTANTIATE_TYPED_TEST_CASE_P(dot_bracket3, alphabet_constexpr, dot_bracket3);

// concepts
TEST(dot_bracket3, concept_check)
{
    EXPECT_TRUE(RnaStructureAlphabet<dot_bracket3>);
    EXPECT_TRUE(RnaStructureAlphabet<dot_bracket3 &>);
    EXPECT_TRUE(RnaStructureAlphabet<dot_bracket3 const>);
    EXPECT_TRUE(RnaStructureAlphabet<dot_bracket3 const &>);

    EXPECT_NE(max_pseudoknot_depth<dot_bracket3>, 0);
}

// assign_char functions
TEST(dot_bracket3, assign_char)
{
    std::vector<char> input
    {
        '.', '(', ')',
        ':', ',', '-', '_', '~', ';',
        '<', '>', '[', ']', '{', '}',
        'H', 'B', 'E', 'G', 'I', 'T', 'S'
    };

    std::vector<dot_bracket3> cmp
    {
        '.'_db3, '('_db3, ')'_db3,
        '.'_db3, '.'_db3, '.'_db3, '.'_db3, '.'_db3, '.'_db3,
        '.'_db3, '.'_db3, '.'_db3, '.'_db3, '.'_db3, '.'_db3,
        '.'_db3, '.'_db3, '.'_db3, '.'_db3, '.'_db3, '.'_db3, '.'_db3
    };

    for (auto [ ch, cm ] : std::view::zip(input, cmp))
        EXPECT_EQ((assign_char_to(ch, dot_bracket3{})), cm);
}

// to_char functions
TEST(dot_bracket3, to_char)
{
    EXPECT_EQ(to_char('.'_db3), '.');
    EXPECT_EQ(to_char('('_db3), '(');
    EXPECT_EQ(to_char(')'_db3), ')');
}

TEST(dot_bracket3, literals)
{
    std::vector<dot_bracket3> vec1;
    vec1.resize(5, '('_db3);
    EXPECT_EQ(vec1, "((((("_db3);

    std::vector<dot_bracket3> vec2{'.'_db3, '('_db3, '('_db3, ')'_db3, ')'_db3, '.'_db3};
    EXPECT_EQ(vec2, ".(())."_db3);
}

TEST(dot_bracket3, rna_structure_properties)
{
    EXPECT_EQ(dot_bracket3::max_pseudoknot_depth, 1);
    EXPECT_TRUE('.'_db3.is_unpaired());
    EXPECT_TRUE('('_db3.is_pair_open());
    EXPECT_TRUE(')'_db3.is_pair_close());

    EXPECT_FALSE('.'_db3.pseudoknot_id());
    EXPECT_EQ('('_db3.pseudoknot_id().value(), 0);
    EXPECT_EQ(')'_db3.pseudoknot_id().value(), 0);
}

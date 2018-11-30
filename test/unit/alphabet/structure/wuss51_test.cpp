// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2018, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <range/v3/view/zip.hpp>

#include <seqan3/alphabet/structure/rna_structure_concept.hpp>
#include <seqan3/alphabet/structure/wuss.hpp>

#include "../alphabet_test_template.hpp"
#include "../alphabet_constexpr_test_template.hpp"

INSTANTIATE_TYPED_TEST_CASE_P(wuss51, alphabet, wuss51);
INSTANTIATE_TYPED_TEST_CASE_P(wuss51, alphabet_constexpr, wuss51);

// assign_char functions
TEST(wuss51, assign_char)
{

    using t = wuss51;
    std::vector<char> input
    {
        '.', '(', ')',
        ':', ',', '-', '_', '~', ';',
        '<', '>', '[', ']', '{', '}',
        'H', 'B', 'E', 'G', 'I', 'T', 'S'
    };

    std::vector<wuss51> cmp
    {
        t::UNPAIRED, t::PAIR_OPEN1, t::PAIR_CLOSE1,
        t::UNPAIRED1, t::UNPAIRED2, t::UNPAIRED3, t::UNPAIRED4, t::UNPAIRED5, t::UNPAIRED6,
        t::PAIR_OPEN, t::PAIR_CLOSE, t::PAIR_OPEN2, t::PAIR_CLOSE2, t::PAIR_OPEN3, t::PAIR_CLOSE3,
        "H"_wuss51.front(),
        "B"_wuss51.front(),
        "E"_wuss51.front(),
        "G"_wuss51.front(),
        "I"_wuss51.front(),
        "T"_wuss51.front(),
        "S"_wuss51.front()
    };

    for (auto [ ch, cm ] : ranges::view::zip(input, cmp))
        EXPECT_EQ((assign_char(wuss51{}, ch)), cm);
}

// to_char functions
TEST(wuss51, to_char)
{
    EXPECT_EQ(to_char(wuss51::UNPAIRED), '.');
    EXPECT_EQ(to_char(wuss51::UNPAIRED1), ':');
    EXPECT_EQ(to_char(wuss51::UNPAIRED2), ',');
    EXPECT_EQ(to_char(wuss51::UNPAIRED3), '-');
    EXPECT_EQ(to_char(wuss51::UNPAIRED4), '_');
    EXPECT_EQ(to_char(wuss51::UNPAIRED5), '~');
    EXPECT_EQ(to_char(wuss51::UNPAIRED6), ';');
    EXPECT_EQ(to_char(wuss51::PAIR_OPEN), '<');
    EXPECT_EQ(to_char(wuss51::PAIR_CLOSE), '>');
    EXPECT_EQ(to_char(wuss51::PAIR_OPEN1), '(');
    EXPECT_EQ(to_char(wuss51::PAIR_CLOSE1), ')');
    EXPECT_EQ(to_char(wuss51::PAIR_OPEN2), '[');
    EXPECT_EQ(to_char(wuss51::PAIR_CLOSE2), ']');
    EXPECT_EQ(to_char(wuss51::PAIR_OPEN3), '{');
    EXPECT_EQ(to_char(wuss51::PAIR_CLOSE3), '}');
}

// concepts
TEST(wuss51, concept_check)
{
    EXPECT_TRUE(rna_structure_concept<wuss51>);
    EXPECT_NE(max_pseudoknot_depth_v<wuss51>, 0);

    EXPECT_TRUE(rna_structure_concept<wuss<>>);  // same as wuss51
    EXPECT_TRUE(rna_structure_concept<wuss<67>>);
}

TEST(wuss51, literals)
{

    std::vector<wuss51> vec1;
    vec1.resize(5, wuss51::PAIR_OPEN);
    EXPECT_EQ(vec1, "<<<<<"_wuss51);

    std::vector<wuss51> vec2{wuss51::UNPAIRED, wuss51::PAIR_OPEN, wuss51::PAIR_OPEN,
                             wuss51::PAIR_CLOSE, wuss51::PAIR_CLOSE, wuss51::UNPAIRED};
    EXPECT_EQ(vec2, ".<<>>."_wuss51);
}

TEST(wuss51, wuss51)
{

    EXPECT_EQ(wuss51::max_pseudoknot_depth, 22);
    std::vector<wuss51> vec = ".:,-_~;<>()[]{}AaBbCcDd"_wuss51;
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

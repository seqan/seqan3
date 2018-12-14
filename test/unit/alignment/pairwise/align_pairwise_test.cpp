// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <type_traits>
#include <utility>

#include <meta/meta.hpp>

#include <seqan3/alignment/configuration/all.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/core/concept/tuple.hpp>
#include <seqan3/range/view/to_char.hpp>

#include <range/v3/view/generate_n.hpp>

using namespace seqan3;

TEST(align_pairwise, single_rng_lvalue)
{

    auto seq1 = "ACGTGATG"_dna4;
    auto seq2 = "AGTGATACT"_dna4;

    auto p = std::make_pair(seq1, seq2);

    {  // the score
        configuration cfg = align_cfg::edit | align_cfg::result{align_cfg::with_score};
        for (auto && res : align_pairwise(p, cfg))
        {
            EXPECT_EQ(res.get_score(), -4);
        }
    }

    {  // the trace
        configuration cfg = align_cfg::edit | align_cfg::result{align_cfg::with_trace};
        for (auto && res : align_pairwise(p, cfg))
        {
            EXPECT_EQ(res.get_score(), -4);
            auto [cmp1, cmp2] = res.get_end_coordinate();
            EXPECT_EQ((std::tie(cmp1, cmp2)), (std::tuple{7, 8}));
            auto && [gap1, gap2] = res.get_alignment();
            EXPECT_EQ(std::string{gap1 | view::to_char}, "ACGTGATG--");
            EXPECT_EQ(std::string{gap2 | view::to_char}, "A-GTGATACT");
        }
    }
}

TEST(align_pairwise, single_view_lvalue)
{

    auto seq1 = "ACGTGATG"_dna4;
    auto seq2 = "AGTGATACT"_dna4;

    // auto p = std::make_pair(seq1, seq2);
    auto v = ranges::view::single(std::tie(seq1, seq2)) | ranges::view::bounded;

    {  // the score
        configuration cfg = align_cfg::edit | align_cfg::result{align_cfg::with_score};
        for (auto && res : align_pairwise(v, cfg))
        {
             EXPECT_EQ(res.get_score(), -4);
        }
    }
    {  // the trace
        configuration cfg = align_cfg::edit | align_cfg::result{align_cfg::with_trace};

        for (auto && res : align_pairwise(v, cfg))
        {
            EXPECT_EQ(res.get_score(), -4);
            auto && [gap1, gap2] = res.get_alignment();
            EXPECT_EQ(std::string{gap1 | view::to_char}, "ACGTGATG--");
            EXPECT_EQ(std::string{gap2 | view::to_char}, "A-GTGATACT");
        }
    }
}

TEST(align_pairwise, multiple_rng_lvalue)
{

    auto seq1 = "ACGTGATG"_dna4;
    auto seq2 = "AGTGATACT"_dna4;

    auto p = std::make_pair(seq1, seq2);
    std::vector<decltype(p)> vec{10, p};

    configuration cfg = align_cfg::edit | align_cfg::result{align_cfg::with_trace};
    for (auto && res : align_pairwise(vec, cfg))
    {
        EXPECT_EQ(res.get_score(), -4);
        auto && [gap1, gap2] = res.get_alignment();
        EXPECT_EQ(std::string{gap1 | view::to_char}, "ACGTGATG--");
        EXPECT_EQ(std::string{gap2 | view::to_char}, "A-GTGATACT");
    }
}

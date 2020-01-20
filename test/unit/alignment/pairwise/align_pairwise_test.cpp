// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
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
#include <seqan3/range/views/to_char.hpp>
#include <seqan3/range/views/to.hpp>
#include <seqan3/std/ranges>

#include <range/v3/view/generate_n.hpp>

using namespace seqan3;

template <typename t>
struct align_pairwise_test : ::testing::Test
{
    // two helper variables to check if the TypeParam contains vectorise.
    using dummy_cfg_t = std::conditional_t<std::is_same_v<void, t>,
                                           decltype(align_cfg::max_error{1}),
                                           t>;
    using config_t = decltype(align_cfg::edit | dummy_cfg_t{});

    static constexpr bool is_vectorised = config_t::template exists<detail::vectorise_tag>();
};

using testing_types = ::testing::Types<void,
                                       align_cfg::parallel>;

TYPED_TEST_SUITE(align_pairwise_test, testing_types, );

template <typename type_param_t, typename seq_t, typename cfg_t>
auto call_alignment(seq_t && seq, cfg_t && cfg)
{
    if constexpr (std::same_as<type_param_t, void>)
    {
        return align_pairwise(std::forward<seq_t>(seq), std::forward<cfg_t>(cfg));
    }
    else
    {
        auto && config = cfg | type_param_t{};
        return align_pairwise(std::forward<seq_t>(seq), std::forward<decltype(config)>(config));
    }
}

TYPED_TEST(align_pairwise_test, single_pair)
{
    auto seq1 = "ACGTGATG"_dna4;
    auto seq2 = "AGTGATACT"_dna4;

    auto p = std::tie(seq1, seq2);

    {  // the score
        configuration cfg = align_cfg::edit | align_cfg::result{with_score};

        for (auto && res : call_alignment<TypeParam>(p, cfg))
        {
            EXPECT_EQ(res.score(), -4.0);
        }
    }

    {  // the alignment
        configuration cfg = align_cfg::edit | align_cfg::result{with_alignment};
        unsigned idx = 0;
        for (auto && res : call_alignment<TypeParam>(p, cfg))
        {
            EXPECT_EQ(res.id(), idx++);
            EXPECT_EQ(res.score(), -4);
            EXPECT_EQ(res.back_coordinate().first, 8u);
            EXPECT_EQ(res.back_coordinate().second, 9u);
            auto && [gap1, gap2] = res.alignment();
            EXPECT_EQ(gap1 | views::to_char | views::to<std::string>, "ACGTGATG--");
            EXPECT_EQ(gap2 | views::to_char | views::to<std::string>, "A-GTGATACT");
        }
    }
}

TYPED_TEST(align_pairwise_test, single_pair_double_score)
{
    if constexpr (!TestFixture::is_vectorised)
    { // not building for vectorised version.
        auto seq1 = "ACGTGATG"_dna4;
        auto seq2 = "AGTGATACT"_dna4;

        auto p = std::tie(seq1, seq2);

        {  // the score
            configuration cfg = align_cfg::edit | align_cfg::result{with_score, using_score_type<double>};

            for (auto && res : call_alignment<TypeParam>(p, cfg))
            {
                EXPECT_EQ(res.score(), -4.0);
                EXPECT_TRUE((std::same_as<decltype(res.score()), double>));
            }
        }

        {  // the alignment
            configuration cfg = align_cfg::edit | align_cfg::result{with_alignment, using_score_type<double>};
            unsigned idx = 0;
            for (auto && res : call_alignment<TypeParam>(p, cfg))
            {
                EXPECT_EQ(res.id(), idx++);
                EXPECT_EQ(res.score(), -4);
                EXPECT_EQ(res.back_coordinate().first, 8u);
                EXPECT_EQ(res.back_coordinate().second, 9u);
                auto && [gap1, gap2] = res.alignment();
                EXPECT_EQ(gap1 | views::to_char | views::to<std::string>, "ACGTGATG--");
                EXPECT_EQ(gap2 | views::to_char | views::to<std::string>, "A-GTGATACT");
            }
        }
    }
}

TYPED_TEST(align_pairwise_test, single_view)
{
    auto seq1 = "ACGTGATG"_dna4;
    auto seq2 = "AGTGATACT"_dna4;

    auto v = std::views::single(std::tie(seq1, seq2)) | std::views::common;

    {  // the score
        configuration cfg = align_cfg::edit | align_cfg::result{with_score};
        for (auto && res : call_alignment<TypeParam>(v, cfg))
        {
             EXPECT_EQ(res.score(), -4);
        }
    }
    {  // the alignment
        configuration cfg = align_cfg::edit | align_cfg::result{with_alignment};

        for (auto && res : call_alignment<TypeParam>(v, cfg))
        {
            EXPECT_EQ(res.score(), -4);
            auto && [gap1, gap2] = res.alignment();
            EXPECT_EQ(gap1 | views::to_char | views::to<std::string>, "ACGTGATG--");
            EXPECT_EQ(gap2 | views::to_char | views::to<std::string>, "A-GTGATACT");
        }
    }
}

TYPED_TEST(align_pairwise_test, collection)
{
    auto seq1 = "ACGTGATG"_dna4;
    auto seq2 = "AGTGATACT"_dna4;

    auto p = std::tie(seq1, seq2);
    std::vector<decltype(p)> vec{10, p};

    configuration cfg = align_cfg::edit | align_cfg::result{with_alignment};
    for (auto && res : call_alignment<TypeParam>(vec, cfg))
    {
        EXPECT_EQ(res.score(), -4);
        auto && [gap1, gap2] = res.alignment();
        EXPECT_EQ(gap1 | views::to_char | views::to<std::string>, "ACGTGATG--");
        EXPECT_EQ(gap2 | views::to_char | views::to<std::string>, "A-GTGATACT");
    }
}

TYPED_TEST(align_pairwise_test, collection_with_double_score_type)
{
    if constexpr (!TestFixture::is_vectorised)
    {
        auto seq1 = "ACGTGATG"_dna4;
        auto seq2 = "AGTGATACT"_dna4;

        auto p = std::tie(seq1, seq2);
        std::vector<decltype(p)> vec{10, p};

        configuration cfg = align_cfg::edit | align_cfg::result{with_alignment, using_score_type<double>};
        for (auto && res : call_alignment<TypeParam>(vec, cfg))
        {
            EXPECT_EQ(res.score(), -4);
            auto && [gap1, gap2] = res.alignment();
            EXPECT_EQ(gap1 | views::to_char | views::to<std::string>, "ACGTGATG--");
            EXPECT_EQ(gap2 | views::to_char | views::to<std::string>, "A-GTGATACT");
        }
    }
}

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

using seqan3::operator""_dna4;

template <typename t>
struct align_pairwise_test : ::testing::Test
{
    // two helper variables to check if the TypeParam contains vectorise.
    using dummy_cfg_t = std::conditional_t<std::is_same_v<void, t>,
                                           decltype(seqan3::align_cfg::max_error{1}),
                                           t>;
    using config_t = decltype(seqan3::align_cfg::edit | dummy_cfg_t{});

    static constexpr bool is_vectorised = config_t::template exists<seqan3::detail::vectorise_tag>();
};

using testing_types = ::testing::Types<void,
                                       seqan3::align_cfg::parallel>;

TYPED_TEST_SUITE(align_pairwise_test, testing_types, );

template <typename type_param_t, typename seq_t, typename cfg_t>
auto call_alignment(seq_t && seq, cfg_t && cfg)
{
    if constexpr (std::same_as<type_param_t, void>)
    {
        return seqan3::align_pairwise(std::forward<seq_t>(seq), std::forward<cfg_t>(cfg));
    }
    else
    {
        auto && config = cfg | type_param_t{};
        return seqan3::align_pairwise(std::forward<seq_t>(seq), std::forward<decltype(config)>(config));
    }
}

TYPED_TEST(align_pairwise_test, single_pair)
{
    auto seq1 = "ACGTGATG"_dna4;
    auto seq2 = "AGTGATACT"_dna4;

    auto p = std::tie(seq1, seq2);

    {  // the score
        seqan3::configuration cfg = seqan3::align_cfg::edit | seqan3::align_cfg::result{seqan3::with_score};

        for (auto && res : call_alignment<TypeParam>(p, cfg))
        {
            EXPECT_EQ(res.score(), -4.0);
        }
    }

    {  // the alignment
        seqan3::configuration cfg = seqan3::align_cfg::edit | seqan3::align_cfg::result{seqan3::with_alignment};
        unsigned idx = 0;
        for (auto && res : call_alignment<TypeParam>(p, cfg))
        {
            EXPECT_EQ(res.id(), idx++);
            EXPECT_EQ(res.score(), -4);
            EXPECT_EQ(res.back_coordinate().first, 8u);
            EXPECT_EQ(res.back_coordinate().second, 9u);
            auto && [gap1, gap2] = res.alignment();
            EXPECT_EQ(gap1 | seqan3::views::to_char | seqan3::views::to<std::string>, "ACGTGATG--");
            EXPECT_EQ(gap2 | seqan3::views::to_char | seqan3::views::to<std::string>, "A-GTGATACT");
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
            seqan3::configuration cfg = seqan3::align_cfg::edit
                                      | seqan3::align_cfg::result{seqan3::with_score, seqan3::using_score_type<double>};

            for (auto && res : call_alignment<TypeParam>(p, cfg))
            {
                EXPECT_EQ(res.score(), -4.0);
                EXPECT_TRUE((std::same_as<decltype(res.score()), double>));
            }
        }

        {  // the alignment
            seqan3::configuration cfg = seqan3::align_cfg::edit |
                                        seqan3::align_cfg::result{seqan3::with_alignment,
                                                                  seqan3::using_score_type<double>};
            unsigned idx = 0;
            for (auto && res : call_alignment<TypeParam>(p, cfg))
            {
                EXPECT_EQ(res.id(), idx++);
                EXPECT_EQ(res.score(), -4);
                EXPECT_EQ(res.back_coordinate().first, 8u);
                EXPECT_EQ(res.back_coordinate().second, 9u);
                auto && [gap1, gap2] = res.alignment();
                EXPECT_EQ(gap1 | seqan3::views::to_char | seqan3::views::to<std::string>, "ACGTGATG--");
                EXPECT_EQ(gap2 | seqan3::views::to_char | seqan3::views::to<std::string>, "A-GTGATACT");
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
        seqan3::configuration cfg = seqan3::align_cfg::edit | seqan3::align_cfg::result{seqan3::with_score};
        for (auto && res : call_alignment<TypeParam>(v, cfg))
        {
             EXPECT_EQ(res.score(), -4);
        }
    }
    {  // the alignment
        seqan3::configuration cfg = seqan3::align_cfg::edit | seqan3::align_cfg::result{seqan3::with_alignment};

        for (auto && res : call_alignment<TypeParam>(v, cfg))
        {
            EXPECT_EQ(res.score(), -4);
            auto && [gap1, gap2] = res.alignment();
            EXPECT_EQ(gap1 | seqan3::views::to_char | seqan3::views::to<std::string>, "ACGTGATG--");
            EXPECT_EQ(gap2 | seqan3::views::to_char | seqan3::views::to<std::string>, "A-GTGATACT");
        }
    }
}

TYPED_TEST(align_pairwise_test, collection)
{
    auto seq1 = "ACGTGATG"_dna4;
    auto seq2 = "AGTGATACT"_dna4;

    auto p = std::tie(seq1, seq2);
    std::vector<decltype(p)> vec{10, p};

    seqan3::configuration cfg = seqan3::align_cfg::edit | seqan3::align_cfg::result{seqan3::with_alignment};
    for (auto && res : call_alignment<TypeParam>(vec, cfg))
    {
        EXPECT_EQ(res.score(), -4);
        auto && [gap1, gap2] = res.alignment();
        EXPECT_EQ(gap1 | seqan3::views::to_char | seqan3::views::to<std::string>, "ACGTGATG--");
        EXPECT_EQ(gap2 | seqan3::views::to_char | seqan3::views::to<std::string>, "A-GTGATACT");
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

        seqan3::configuration cfg = seqan3::align_cfg::edit |
                                    seqan3::align_cfg::result{seqan3::with_alignment, seqan3::using_score_type<double>};
        for (auto && res : call_alignment<TypeParam>(vec, cfg))
        {
            EXPECT_EQ(res.score(), -4);
            auto && [gap1, gap2] = res.alignment();
            EXPECT_EQ(gap1 | seqan3::views::to_char | seqan3::views::to<std::string>, "ACGTGATG--");
            EXPECT_EQ(gap2 | seqan3::views::to_char | seqan3::views::to<std::string>, "A-GTGATACT");
        }
    }
}

TYPED_TEST(align_pairwise_test, bug_1598)
{
    // https://github.com/seqan/seqan3/issues/1598
    using seqan3::operator""_dna4;
    auto s1 = "TTACGTACGGACTAGCTACAACATTACGGACTAC"_dna4;
    auto g  = "GGACGACATGACGTACGACTTTACGTACGACTAGC"_dna4;
    auto s2 = g | std::views::drop(2);

    // Configure the alignment kernel.
    seqan3::configuration cfg = seqan3::align_cfg::mode{seqan3::global_alignment} |
                                seqan3::align_cfg::scoring{seqan3::nucleotide_scoring_scheme{}} |
                                seqan3::align_cfg::result{seqan3::with_alignment};

    // Invoke the pairwise alignment which returns a lazy range over alignment results.
    auto results = seqan3::align_pairwise(std::tie(s1, s2), cfg);
}

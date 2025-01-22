// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <ranges>
#include <stdexcept>
#include <string_view>
#include <type_traits>
#include <utility>

#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/views/to_char.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/test/expect_same_type.hpp>
#include <seqan3/utility/tuple/concept.hpp>

using seqan3::operator""_dna4;

template <typename t>
struct align_pairwise_test : ::testing::Test
{
    // two helper variables to check if the TypeParam contains vectorised.
    using dummy_cfg_t = std::conditional_t<std::is_same_v<void, t>, decltype(seqan3::align_cfg::min_score{-1}), t>;
    using config_t = decltype(seqan3::align_cfg::method_global{} | seqan3::align_cfg::edit_scheme | dummy_cfg_t{});

    static constexpr bool is_vectorised = config_t::template exists<seqan3::align_cfg::vectorised>();
};

using testing_types = ::testing::Types<void, seqan3::align_cfg::parallel>;

TYPED_TEST_SUITE(align_pairwise_test, testing_types, );

template <typename type_param_t, typename seq_t, typename cfg_t>
auto call_alignment(seq_t && seq, cfg_t && cfg)
{
    if constexpr (std::same_as<type_param_t, void>)
    {
        return seqan3::align_pairwise(std::forward<seq_t>(seq), std::forward<cfg_t>(cfg));
    }
    else if constexpr (std::same_as<type_param_t, seqan3::align_cfg::parallel>)
    {
        auto && config = cfg | seqan3::align_cfg::parallel{4};
        return seqan3::align_pairwise(std::forward<seq_t>(seq), std::forward<decltype(config)>(config));
    }
}

TYPED_TEST(align_pairwise_test, single_pair)
{
    using namespace std::literals;

    auto seq1 = "ACGTGATG"_dna4;
    auto seq2 = "AGTGATACT"_dna4;

    // use std::tie
    auto p1 = std::tie(seq1, seq2);

    { // the score
        seqan3::configuration cfg =
            seqan3::align_cfg::method_global{} | seqan3::align_cfg::edit_scheme | seqan3::align_cfg::output_score{};

        for (auto && res : call_alignment<TypeParam>(p1, cfg))
        {
            EXPECT_EQ(res.score(), -4.0);
        }
    }

    { // with everything (default if no output is specified)
        seqan3::configuration cfg = seqan3::align_cfg::method_global{} | seqan3::align_cfg::edit_scheme;
        unsigned idx = 0;
        for (auto && res : call_alignment<TypeParam>(p1, cfg))
        {
            EXPECT_EQ(res.sequence1_id(), idx);
            EXPECT_EQ(res.sequence2_id(), idx++);
            EXPECT_EQ(res.score(), -4);
            EXPECT_EQ(res.sequence1_end_position(), 8u);
            EXPECT_EQ(res.sequence2_end_position(), 9u);
            auto && [gap1, gap2] = res.alignment();
            EXPECT_RANGE_EQ(gap1 | seqan3::views::to_char, "ACGTGATG--"sv);
            EXPECT_RANGE_EQ(gap2 | seqan3::views::to_char, "A-GTGATACT"sv);
        }
    }

    // use std::pair
    std::pair sequences{"ACGTGATG"_dna4, "AGTGATACT"_dna4};

    { // the score
        seqan3::configuration cfg =
            seqan3::align_cfg::method_global{} | seqan3::align_cfg::edit_scheme | seqan3::align_cfg::output_score{};

        for (auto && res : call_alignment<TypeParam>(sequences, cfg))
        {
            EXPECT_EQ(res.score(), -4.0);
        }
    }

    { // with everything (default if no output is specified)
        seqan3::configuration cfg = seqan3::align_cfg::method_global{} | seqan3::align_cfg::edit_scheme;
        unsigned idx = 0;
        for (auto && res : call_alignment<TypeParam>(sequences, cfg))
        {
            EXPECT_EQ(res.sequence1_id(), idx);
            EXPECT_EQ(res.sequence2_id(), idx++);
            EXPECT_EQ(res.score(), -4);
            EXPECT_EQ(res.sequence1_end_position(), 8u);
            EXPECT_EQ(res.sequence2_end_position(), 9u);
            auto && [gap1, gap2] = res.alignment();
            EXPECT_RANGE_EQ(gap1 | seqan3::views::to_char, "ACGTGATG--"sv);
            EXPECT_RANGE_EQ(gap2 | seqan3::views::to_char, "A-GTGATACT"sv);
        }
    }

    // use make_pair
    auto p2 = std::make_pair(seq1, seq2);

    { // the score
        seqan3::configuration cfg =
            seqan3::align_cfg::method_global{} | seqan3::align_cfg::edit_scheme | seqan3::align_cfg::output_score{};

        for (auto && res : call_alignment<TypeParam>(p2, cfg))
        {
            EXPECT_EQ(res.score(), -4.0);
        }
    }

    { // with everything (default if no output is specified)
        seqan3::configuration cfg = seqan3::align_cfg::method_global{} | seqan3::align_cfg::edit_scheme;
        unsigned idx = 0;
        for (auto && res : call_alignment<TypeParam>(p2, cfg))
        {
            EXPECT_EQ(res.sequence1_id(), idx);
            EXPECT_EQ(res.sequence2_id(), idx++);
            EXPECT_EQ(res.score(), -4);
            EXPECT_EQ(res.sequence1_end_position(), 8u);
            EXPECT_EQ(res.sequence2_end_position(), 9u);
            auto && [gap1, gap2] = res.alignment();
            EXPECT_RANGE_EQ(gap1 | seqan3::views::to_char, "ACGTGATG--"sv);
            EXPECT_RANGE_EQ(gap2 | seqan3::views::to_char, "A-GTGATACT"sv);
        }
    }
}

TYPED_TEST(align_pairwise_test, single_pair_double_score)
{
    using namespace std::literals;

    if constexpr (!TestFixture::is_vectorised)
    { // not building for vectorised version.
        auto seq1 = "ACGTGATG"_dna4;
        auto seq2 = "AGTGATACT"_dna4;

        auto p = std::tie(seq1, seq2);

        { // the score
            seqan3::configuration cfg = seqan3::align_cfg::method_global{} | seqan3::align_cfg::edit_scheme
                                      | seqan3::align_cfg::output_score{} | seqan3::align_cfg::score_type<double>{};

            for (auto && res : call_alignment<TypeParam>(p, cfg))
            {
                EXPECT_EQ(res.score(), -4.0);
                EXPECT_SAME_TYPE(decltype(res.score()), double);
            }
        }

        { // with everything (default if no output is specified)
            seqan3::configuration cfg = seqan3::align_cfg::method_global{} | seqan3::align_cfg::edit_scheme
                                      | seqan3::align_cfg::output_alignment{} | seqan3::align_cfg::output_end_position{}
                                      | seqan3::align_cfg::output_sequence1_id{}
                                      | seqan3::align_cfg::output_sequence2_id{} | seqan3::align_cfg::output_score{}
                                      | seqan3::align_cfg::score_type<double>{};

            for (auto && res : call_alignment<TypeParam>(p, cfg))
            {
                EXPECT_EQ(res.score(), -4);
                EXPECT_EQ(res.sequence1_end_position(), 8u);
                EXPECT_EQ(res.sequence2_end_position(), 9u);
                auto && [gap1, gap2] = res.alignment();
                EXPECT_RANGE_EQ(gap1 | seqan3::views::to_char, "ACGTGATG--"sv);
                EXPECT_RANGE_EQ(gap2 | seqan3::views::to_char, "A-GTGATACT"sv);
            }
        }
    }
}

TYPED_TEST(align_pairwise_test, single_view)
{
    using namespace std::literals;

    auto seq1 = "ACGTGATG"_dna4;
    auto seq2 = "AGTGATACT"_dna4;

    auto v = std::views::single(std::tie(seq1, seq2)) | std::views::common;

    { // the score
        seqan3::configuration cfg =
            seqan3::align_cfg::method_global{} | seqan3::align_cfg::edit_scheme | seqan3::align_cfg::output_score{};
        for (auto && res : call_alignment<TypeParam>(v, cfg))
        {
            EXPECT_EQ(res.score(), -4);
        }
    }

    { // with everything (default if no output is specified)
        seqan3::configuration cfg = seqan3::align_cfg::method_global{} | seqan3::align_cfg::edit_scheme;

        for (auto && res : call_alignment<TypeParam>(v, cfg))
        {
            EXPECT_EQ(res.score(), -4);
            auto && [gap1, gap2] = res.alignment();
            EXPECT_RANGE_EQ(gap1 | seqan3::views::to_char, "ACGTGATG--"sv);
            EXPECT_RANGE_EQ(gap2 | seqan3::views::to_char, "A-GTGATACT"sv);
        }
    }
}

TYPED_TEST(align_pairwise_test, collection)
{
    using namespace std::literals;

    auto seq1 = "ACGTGATG"_dna4;
    auto seq2 = "AGTGATACT"_dna4;

    auto p = std::tie(seq1, seq2);
    std::vector<decltype(p)> vec{10, p};

    seqan3::configuration cfg = seqan3::align_cfg::method_global{} | seqan3::align_cfg::edit_scheme
                              | seqan3::align_cfg::output_score{} | seqan3::align_cfg::output_alignment{};
    for (auto && res : call_alignment<TypeParam>(vec, cfg))
    {
        EXPECT_EQ(res.score(), -4);
        auto && [gap1, gap2] = res.alignment();
        EXPECT_RANGE_EQ(gap1 | seqan3::views::to_char, "ACGTGATG--"sv);
        EXPECT_RANGE_EQ(gap2 | seqan3::views::to_char, "A-GTGATACT"sv);
    }
}

TYPED_TEST(align_pairwise_test, collection_with_double_score_type)
{
    using namespace std::literals;

    if constexpr (!TestFixture::is_vectorised)
    {
        auto seq1 = "ACGTGATG"_dna4;
        auto seq2 = "AGTGATACT"_dna4;

        auto p = std::tie(seq1, seq2);
        std::vector<decltype(p)> vec{10, p};

        seqan3::configuration cfg = seqan3::align_cfg::method_global{} | seqan3::align_cfg::edit_scheme
                                  | seqan3::align_cfg::output_alignment{} | seqan3::align_cfg::output_score{}
                                  | seqan3::align_cfg::score_type<double>{};

        for (auto && res : call_alignment<TypeParam>(vec, cfg))
        {
            EXPECT_EQ(res.score(), -4);
            auto && [gap1, gap2] = res.alignment();
            EXPECT_RANGE_EQ(gap1 | seqan3::views::to_char, "ACGTGATG--"sv);
            EXPECT_RANGE_EQ(gap2 | seqan3::views::to_char, "A-GTGATACT"sv);
        }
    }
}

TYPED_TEST(align_pairwise_test, bug_1598)
{
    // https://github.com/seqan/seqan3/issues/1598
    using seqan3::operator""_dna4;
    auto s1 = "TTACGTACGGACTAGCTACAACATTACGGACTAC"_dna4;
    auto g = "GGACGACATGACGTACGACTTTACGTACGACTAGC"_dna4;
    auto s2 = g | std::views::drop(2);

    // Configure the alignment kernel.
    seqan3::configuration cfg = seqan3::align_cfg::method_global{}
                              | seqan3::align_cfg::scoring_scheme{seqan3::nucleotide_scoring_scheme{}}
                              | seqan3::align_cfg::output_alignment{};

    // Invoke the pairwise alignment which returns a lazy range over alignment results.
    auto results = seqan3::align_pairwise(std::tie(s1, s2), cfg);
}

TEST(align_pairwise_test, parallel_without_parameter)
{
    auto seq1 = "ACGTGATG"_dna4;
    auto seq2 = "AGTGATACT"_dna4;
    seqan3::configuration cfg = seqan3::align_cfg::method_global{} | seqan3::align_cfg::edit_scheme
                              | seqan3::align_cfg::output_score{} | seqan3::align_cfg::parallel{};

    EXPECT_THROW(seqan3::align_pairwise(std::tie(seq1, seq2), cfg), std::runtime_error);
}

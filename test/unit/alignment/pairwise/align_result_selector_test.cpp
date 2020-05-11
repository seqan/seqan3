// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <list>
#include <type_traits>
#include <utility>
#include <vector>

#include <meta/meta.hpp>

#include <seqan3/alignment/configuration/all.hpp>
#include <seqan3/alignment/configuration/align_config_result.hpp>
#include <seqan3/alignment/pairwise/align_result_selector.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/concept/tuple.hpp>
#include <seqan3/range/views/type_reduce.hpp>

TEST(alignment_selector, align_result_selector_with_list)
{
    using seq1_t = std::vector<seqan3::dna4>;
    using seq2_t = std::list<seqan3::dna4>;

    {
        auto cfg = seqan3::align_cfg::edit | seqan3::align_cfg::result{seqan3::with_score};
        using _t = seqan3::alignment_result<typename seqan3::detail::align_result_selector<seq1_t,
                                                                                           seq2_t,
                                                                                           decltype(cfg)>::type>;

        EXPECT_TRUE((std::is_same_v<decltype(std::declval<_t>().id()), uint32_t>));
        EXPECT_TRUE((std::is_same_v<decltype(std::declval<_t>().score()), int32_t>));
    }

    {
        auto cfg = seqan3::align_cfg::edit | seqan3::align_cfg::result{seqan3::with_alignment};
        using _t = seqan3::alignment_result<typename seqan3::detail::align_result_selector<seq1_t,
                                                                                           seq2_t,
                                                                                           decltype(cfg)>::type>;

        using gapped_seq1_t = seqan3::gap_decorator<seqan3::type_reduce_view<std::vector<seqan3::dna4> &>>;
        using gapped_seq2_t = std::vector<seqan3::gapped<seqan3::dna4>>;

        EXPECT_TRUE((std::is_same_v<decltype(std::declval<_t>().id()), uint32_t>));
        EXPECT_TRUE((std::is_same_v<decltype(std::declval<_t>().score()), int32_t>));
        EXPECT_TRUE((std::is_same_v<decltype(std::declval<_t>().back_coordinate()),
                                    seqan3::alignment_coordinate const &>));
        EXPECT_TRUE((std::is_same_v<decltype(std::declval<_t>().front_coordinate()),
                                    seqan3::alignment_coordinate const &>));
        EXPECT_TRUE((std::is_same_v<decltype(std::declval<_t>().alignment()),
                                    std::tuple<gapped_seq1_t, gapped_seq2_t> const &>));
    }
}

TEST(alignment_selector, align_result_selector_using_score_type)
{
    using seq1_t = std::vector<seqan3::dna4>;
    using seq2_t = std::list<seqan3::dna4>;

    auto cfg = seqan3::align_cfg::edit | seqan3::align_cfg::result{seqan3::with_back_coordinate,
                                                                   seqan3::using_score_type<double>};
    using _t = seqan3::alignment_result<typename seqan3::detail::align_result_selector<seq1_t,
                                                                                       seq2_t,
                                                                                       decltype(cfg)>::type>;

    EXPECT_TRUE((std::is_same_v<decltype(std::declval<_t>().score()), double>));
    EXPECT_TRUE((std::is_same_v<decltype(std::declval<_t>().back_coordinate()), seqan3::alignment_coordinate const &>));
}

TEST(alignment_selector, align_result_selector_with_vector)
{
    using seq_t = std::vector<seqan3::dna4>;

    auto cfg = seqan3::align_cfg::edit | seqan3::align_cfg::result{seqan3::with_alignment};
    using _t = seqan3::alignment_result<typename seqan3::detail::align_result_selector<seq_t,
                                                                                       seq_t,
                                                                                       decltype(cfg)>::type>;

    using gapped_seq_t = seqan3::gap_decorator<seqan3::type_reduce_view<std::vector<seqan3::dna4> &>>;

    EXPECT_TRUE((std::is_same_v<decltype(std::declval<_t>().id()), uint32_t>));
    EXPECT_TRUE((std::is_same_v<decltype(std::declval<_t>().score()), int32_t>));
    EXPECT_TRUE((std::is_same_v<decltype(std::declval<_t>().back_coordinate()),
                                seqan3::alignment_coordinate const &>));
    EXPECT_TRUE((std::is_same_v<decltype(std::declval<_t>().front_coordinate()),
                                seqan3::alignment_coordinate const &>));
    EXPECT_TRUE((std::is_same_v<decltype(std::declval<_t>().alignment()),
                                std::tuple<gapped_seq_t, gapped_seq_t> const &>));
}

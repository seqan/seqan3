// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <list>
#include <type_traits>
#include <utility>
#include <vector>

#include <meta/meta.hpp>

#include <seqan3/alignment/configuration/align_config_edit.hpp>
#include <seqan3/alignment/configuration/align_config_result.hpp>
#include <seqan3/alignment/pairwise/alignment_selector.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/core/concept/tuple.hpp>
#include <seqan3/range/view/persist.hpp>
#include <seqan3/range/view/to_char.hpp>

using namespace seqan3;

TEST(alignment_selector, determine_result_type)
{
    using seq1_t = std::vector<dna4>;
    using seq2_t = std::list<dna4>;

    { // test case I
        auto cfg = align_cfg::edit;
        using _t = typename detail::determine_result_type<seq1_t, seq2_t, decltype(cfg)>::type;

        EXPECT_TRUE((std::is_same_v<decltype(std::declval<_t>().get_id()), uint32_t>));
        EXPECT_TRUE((std::is_same_v<decltype(std::declval<_t>().get_score()), int32_t>));
    }

    { // test case II
        auto cfg = align_cfg::edit | align_cfg::result{align_cfg::with_score};
        using _t = typename detail::determine_result_type<seq1_t, seq2_t, decltype(cfg)>::type;

        EXPECT_TRUE((std::is_same_v<decltype(std::declval<_t>().get_id()), uint32_t>));
        EXPECT_TRUE((std::is_same_v<decltype(std::declval<_t>().get_score()), int32_t>));
    }

    { // test case III
        auto cfg = align_cfg::edit | align_cfg::result{align_cfg::with_trace};
        using _t = typename detail::determine_result_type<seq1_t, seq2_t, decltype(cfg)>::type;

        using gapped_seq1_t = std::vector<gapped<dna4>>;
        using gapped_seq2_t = std::vector<gapped<dna4>>;

        EXPECT_TRUE((std::is_same_v<decltype(std::declval<_t>().get_id()), uint32_t>));
        EXPECT_TRUE((std::is_same_v<decltype(std::declval<_t>().get_score()), int32_t>));
        EXPECT_TRUE((std::is_same_v<decltype(std::declval<_t>().get_end_coordinate()),
                                    detail::alignment_coordinate const &>));
        EXPECT_TRUE((std::is_same_v<decltype(std::declval<_t>().get_begin_coordinate()),
                                    detail::alignment_coordinate const &>));
        EXPECT_TRUE((std::is_same_v<decltype(std::declval<_t>().get_alignment()),
                                    std::pair<gapped_seq1_t, gapped_seq2_t> const &>));
    }
}

TEST(alignment_selector, select)
{
    using seq1_t = std::vector<dna4>;
    using seq2_t = std::vector<dna4>;
    using seq_t = std::pair<seq1_t, seq2_t>;

    { // test case I
        auto cfg = align_cfg::edit;
        using _t = detail::alignment_selector<seq_t, decltype(cfg)>;

        using res_t = typename detail::determine_result_type<seq1_t, seq2_t, decltype(cfg)>::type;

        seq_t s{};
        EXPECT_TRUE((std::is_same_v<res_t, typename _t::result_type>));
        EXPECT_TRUE((std::is_same_v<decltype(_t{cfg}.select(s)),
                                    std::function<res_t(res_t &)>>));
    }
}

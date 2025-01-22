// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <list>
#include <type_traits>
#include <utility>
#include <vector>

#include <seqan3/alignment/configuration/align_config_edit.hpp>
#include <seqan3/alignment/configuration/align_config_method.hpp>
#include <seqan3/alignment/configuration/align_config_output.hpp>
#include <seqan3/alignment/configuration/align_config_score_type.hpp>
#include <seqan3/alignment/pairwise/align_result_selector.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/utility/tuple/concept.hpp>
#include <seqan3/utility/views/type_reduce.hpp>

struct alignment_selector_test : public ::testing::Test
{
    using seq1_t = std::vector<seqan3::dna4>;
    using seq2_t = std::list<seqan3::dna4>;

    using gapped_seq1_t = seqan3::gap_decorator<seqan3::type_reduce_t<std::vector<seqan3::dna4> &>>;
    using gapped_seq2_t = std::vector<seqan3::gapped<seqan3::dna4>>;

    using alignment_t = std::tuple<gapped_seq1_t, gapped_seq2_t>;

    template <typename config_t>
    using alignment_result_value_t = typename seqan3::detail::align_result_selector<seq1_t, seq2_t, config_t>::type;

    template <typename config_t>
    using alignment_result_t = seqan3::alignment_result<alignment_result_value_t<config_t>>;

    static constexpr seqan3::configuration base_config = []()
    {
        return seqan3::align_cfg::method_global{} | seqan3::align_cfg::edit_scheme;
    }();
};

TEST_F(alignment_selector_test, align_result_selector_all)
{
    auto cfg = base_config | seqan3::align_cfg::output_score{} | seqan3::align_cfg::output_begin_position{}
             | seqan3::align_cfg::output_end_position{} | seqan3::align_cfg::output_alignment{}
             | seqan3::align_cfg::output_sequence1_id{} | seqan3::align_cfg::output_sequence2_id{};

    using result_t = alignment_result_t<decltype(cfg)>;

    EXPECT_TRUE((std::is_same_v<decltype(std::declval<result_t>().sequence1_id()), uint32_t>));
    EXPECT_TRUE((std::is_same_v<decltype(std::declval<result_t>().sequence2_id()), uint32_t>));
    EXPECT_TRUE((std::is_same_v<decltype(std::declval<result_t>().score()), int32_t>));
    EXPECT_TRUE((std::is_same_v<decltype(std::declval<result_t>().sequence1_end_position()), size_t>));
    EXPECT_TRUE((std::is_same_v<decltype(std::declval<result_t>().sequence2_end_position()), size_t>));
    EXPECT_TRUE((std::is_same_v<decltype(std::declval<result_t>().sequence1_begin_position()), size_t>));
    EXPECT_TRUE((std::is_same_v<decltype(std::declval<result_t>().sequence2_begin_position()), size_t>));
    EXPECT_TRUE((std::is_same_v<decltype(std::declval<result_t>().alignment()), alignment_t const &>));
}

TEST_F(alignment_selector_test, align_result_selector_using_score_type)
{
    auto cfg = seqan3::align_cfg::method_global{} | // TODO: Change to base_config once using_score_type has been fixed.
               seqan3::align_cfg::edit_scheme | seqan3::align_cfg::output_score{}
             | seqan3::align_cfg::output_end_position{} | seqan3::align_cfg::score_type<double>{};

    using result_t = alignment_result_t<decltype(cfg)>;

    EXPECT_TRUE((std::is_same_v<decltype(std::declval<result_t>().score()), double>));
    EXPECT_TRUE((std::is_same_v<decltype(std::declval<result_t>().sequence1_end_position()), size_t>));
    EXPECT_TRUE((std::is_same_v<decltype(std::declval<result_t>().sequence2_end_position()), size_t>));
}

TEST_F(alignment_selector_test, output_score_only)
{
    auto cfg = base_config | seqan3::align_cfg::output_score{};

    using result_t = alignment_result_t<decltype(cfg)>;

    EXPECT_TRUE((std::is_same_v<decltype(std::declval<result_t>().score()), int32_t>));
}

TEST_F(alignment_selector_test, output_end_positions_only)
{
    auto cfg = base_config | seqan3::align_cfg::output_end_position{};

    using result_t = alignment_result_t<decltype(cfg)>;

    EXPECT_TRUE((std::is_same_v<decltype(std::declval<result_t>().sequence1_end_position()), size_t>));
    EXPECT_TRUE((std::is_same_v<decltype(std::declval<result_t>().sequence2_end_position()), size_t>));
}

TEST_F(alignment_selector_test, output_begin_positions_only)
{
    auto cfg = base_config | seqan3::align_cfg::output_begin_position{};

    using result_t = alignment_result_t<decltype(cfg)>;

    EXPECT_TRUE((std::is_same_v<decltype(std::declval<result_t>().sequence1_begin_position()), size_t>));
    EXPECT_TRUE((std::is_same_v<decltype(std::declval<result_t>().sequence2_begin_position()), size_t>));
}

TEST_F(alignment_selector_test, output_alignment_only)
{
    auto cfg = base_config | seqan3::align_cfg::output_alignment{};

    using result_t = alignment_result_t<decltype(cfg)>;

    EXPECT_TRUE((std::is_same_v<decltype(std::declval<result_t>().alignment()), alignment_t const &>));
}

TEST_F(alignment_selector_test, output_sequence1_id_only)
{
    auto cfg = base_config | seqan3::align_cfg::output_sequence1_id{};

    using result_t = alignment_result_t<decltype(cfg)>;

    EXPECT_TRUE((std::is_same_v<decltype(std::declval<result_t>().sequence1_id()), uint32_t>));
}

TEST_F(alignment_selector_test, output_sequence2_id_only)
{
    auto cfg = base_config | seqan3::align_cfg::output_sequence2_id{};

    using result_t = alignment_result_t<decltype(cfg)>;

    EXPECT_TRUE((std::is_same_v<decltype(std::declval<result_t>().sequence2_id()), uint32_t>));
}

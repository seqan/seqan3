// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <functional>
#include <type_traits>

#include <seqan3/alignment/configuration/align_config_sequence_ends.hpp>

TEST(align_config_sequence_ends, constructor)
{
    EXPECT_TRUE((std::is_default_constructible_v<seqan3::detail::align_config_sequence_ends<>>));
}

TEST(align_config_sequence_ends, on_align_config)
{
    struct bar
    {
        int value;
    };

    using sequence_ends_align_cfg = seqan3::detail::on_align_config<seqan3::align_cfg::id::sequence_ends>;

    EXPECT_FALSE((sequence_ends_align_cfg::invoke<bar>::value));
    EXPECT_TRUE((sequence_ends_align_cfg::invoke<seqan3::detail::align_config_sequence_ends<>>::value));
}

TEST(align_config_sequence_ends, align_config_type_to_id)
{
    EXPECT_EQ(seqan3::detail::align_config_type_to_id_v<seqan3::detail::align_config_sequence_ends<>>,
              seqan3::align_cfg::id::sequence_ends);
    EXPECT_EQ(seqan3::detail::align_config_type_to_id<seqan3::detail::align_config_sequence_ends<>>::value,
              seqan3::align_cfg::id::sequence_ends);
}

TEST(align_config_sequence_ends, invoke)
{
    auto cfg = std::invoke(seqan3::align_cfg::sequence_ends<>(free_ends_at::seq1), seqan3::detail::configuration<>{});

    EXPECT_NE(free_ends_at::none, get<0>(cfg).value & free_ends_at::seq1_front);
    EXPECT_EQ(free_ends_at::none, get<0>(cfg).value & free_ends_at::seq2_front);
    EXPECT_NE(free_ends_at::none, get<0>(cfg).value & free_ends_at::seq1_back);
    EXPECT_EQ(free_ends_at::none, get<0>(cfg).value & free_ends_at::seq2_back);

    EXPECT_TRUE((std::is_same_v<remove_cvref_t<decltype(cfg)>,
                                seqan3::detail::configuration<seqan3::detail::align_config_sequence_ends_deferred>>));
}

TEST(align_config_sequence_ends, invoke_static)
{
    auto cfg = std::invoke(seqan3::align_cfg::sequence_ends<free_ends_at::seq1>(), seqan3::detail::configuration<>{});

    EXPECT_NE(free_ends_at::none, get<0>(cfg).value & free_ends_at::seq1_front);
    EXPECT_EQ(free_ends_at::none, get<0>(cfg).value & free_ends_at::seq2_front);
    EXPECT_NE(free_ends_at::none, get<0>(cfg).value & free_ends_at::seq1_back);
    EXPECT_EQ(free_ends_at::none, get<0>(cfg).value & free_ends_at::seq2_back);

    using free_ends_align_config = seqan3::detail::align_config_sequence_ends<free_ends_at::seq1>>;

    EXPECT_TRUE((std::is_same_v<remove_cvref_t<decltype(cfg)>, seqan3::detail::configuration<free_ends_align_config>));
}

TEST(align_config_sequence_ends, get_by_enum)
{
    {
        seqan3::detail::configuration cfg = seqan3::align_cfg::sequence_ends<>(free_ends_at::seq1_back |
                                                                               free_ends_at::seq2_front);

        EXPECT_EQ(free_ends_at::none, get<0>(cfg).value & free_ends_at::seq1_front);
        EXPECT_NE(free_ends_at::none, get<0>(cfg).value & free_ends_at::seq2_front);
        EXPECT_NE(free_ends_at::none, get<0>(cfg).value & free_ends_at::seq1_back);
        EXPECT_EQ(free_ends_at::none, get<0>(cfg).value & free_ends_at::seq2_back);

        EXPECT_TRUE((std::is_same_v<decltype(std::get<seqan3::align_cfg::id::sequence_ends>(cfg)),
                                    free_ends_at &>));
    }

    {
        seqan3::detail::configuration<seqan3::detail::align_config_sequence_ends_deferred> const c_cfg =
            seqan3::detail::configuration{seqan3::align_cfg::sequence_ends<>(free_ends_at::seq1_back |
                                                                             free_ends_at::seq2_front)};

        EXPECT_EQ(free_ends_at::none, get<0>(c_cfg).value & free_ends_at::seq1_front);
        EXPECT_NE(free_ends_at::none, get<0>(c_cfg).value & free_ends_at::seq2_front);
        EXPECT_NE(free_ends_at::none, get<0>(c_cfg).value & free_ends_at::seq1_back);
        EXPECT_EQ(free_ends_at::none, get<0>(c_cfg).value & free_ends_at::seq2_back);

        EXPECT_TRUE((std::is_same_v<decltype(std::get<seqan3::align_cfg::id::sequence_ends>(c_cfg)),
                                    free_ends_at const &>));
    }

    {
        seqan3::detail::configuration cfg = seqan3::align_cfg::sequence_ends<>(free_ends_at::seq1_back |
                                                                               free_ends_at::seq2_front);

        EXPECT_EQ(free_ends_at::none, get<0>(cfg).value & free_ends_at::seq1_front);
        EXPECT_NE(free_ends_at::none, get<0>(cfg).value & free_ends_at::seq2_front);
        EXPECT_NE(free_ends_at::none, get<0>(cfg).value & free_ends_at::seq1_back);
        EXPECT_EQ(free_ends_at::none, get<0>(cfg).value & free_ends_at::seq2_back);

        EXPECT_TRUE((std::is_same_v<decltype(std::get<seqan3::align_cfg::id::sequence_ends>(std::move(cfg))),
                                    free_ends_at &&>));
    }

    {
        seqan3::detail::configuration<seqan3::detail::align_config_sequence_ends_deferred> const c_cfg =
            seqan3::detail::configuration{seqan3::align_cfg::sequence_ends<>(free_ends_at::seq1_back |
                                                                             free_ends_at::seq2_front)};

        EXPECT_EQ(free_ends_at::none, get<0>(c_cfg).value & free_ends_at::seq1_front);
        EXPECT_NE(free_ends_at::none, get<0>(c_cfg).value & free_ends_at::seq2_front);
        EXPECT_NE(free_ends_at::none, get<0>(c_cfg).value & free_ends_at::seq1_back);
        EXPECT_EQ(free_ends_at::none, get<0>(c_cfg).value & free_ends_at::seq2_back);

        EXPECT_TRUE((std::is_same_v<decltype(std::get<seqan3::align_cfg::id::sequence_ends>(std::move(c_cfg))),
                                    free_ends_at const &&>));
    }
}

TEST(align_config_sequence_ends, free_ends_enum_all_and_none)
{
    seqan3::seqan3::detail::align_config_sequence_ends<free_ends_at::all> cfg_all;
    EXPECT_NE(free_ends_at::none, cfg_all.value & free_ends_at::seq1_front);
    EXPECT_NE(free_ends_at::none, cfg_all.value & free_ends_at::seq2_front);
    EXPECT_NE(free_ends_at::none, cfg_all.value & free_ends_at::seq1_back);
    EXPECT_NE(free_ends_at::none, cfg_all.value & free_ends_at::seq2_back);

    seqan3::seqan3::detail::align_config_sequence_ends cfg_none;
    EXPECT_EQ(free_ends_at::none, cfg_none.value & free_ends_at::seq1_front);
    EXPECT_EQ(free_ends_at::none, cfg_none.value & free_ends_at::seq2_front);
    EXPECT_EQ(free_ends_at::none, cfg_none.value & free_ends_at::seq1_back);
    EXPECT_EQ(free_ends_at::none, cfg_none.value & free_ends_at::seq2_back);
}

TEST(align_config_sequence_ends, invoke_deferred)
{
    seqan3::detail::configuration cfg = seqan3::align_cfg::sequence_ends<>(free_ends_at::seq1);

    auto call_on_site = [] (auto && new_cfg)
    {
        return std::get<0>(new_cfg).value;
    };

    EXPECT_EQ((std::invoke(std::get<0>(cfg), call_on_site, cfg)), free_ends_at::seq1);
}

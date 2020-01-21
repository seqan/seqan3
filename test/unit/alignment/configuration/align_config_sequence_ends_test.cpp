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

using namespace seqan3;

TEST(align_config_sequence_ends, constructor)
{
    EXPECT_TRUE((std::is_default_constructible_v<detail::align_config_sequence_ends<>>));
}

TEST(align_config_sequence_ends, on_align_config)
{
    struct bar
    {
        int value;
    };

    EXPECT_FALSE((detail::on_align_config<align_cfg::id::sequence_ends>::invoke<bar>::value));
    EXPECT_TRUE((detail::on_align_config<align_cfg::id::sequence_ends>::invoke
        <detail::align_config_sequence_ends<>>::value));
}

TEST(align_config_sequence_ends, align_config_type_to_id)
{
    EXPECT_EQ(detail::align_config_type_to_id_v<detail::align_config_sequence_ends<>>, align_cfg::id::sequence_ends);
    EXPECT_EQ(detail::align_config_type_to_id<detail::align_config_sequence_ends<>>::value,
              align_cfg::id::sequence_ends);
}

TEST(align_config_sequence_ends, invoke)
{
    auto cfg = std::invoke(align_cfg::sequence_ends<>(free_ends_at::seq1), detail::configuration<>{});

    EXPECT_NE(free_ends_at::none, get<0>(cfg).value & free_ends_at::seq1_front);
    EXPECT_EQ(free_ends_at::none, get<0>(cfg).value & free_ends_at::seq2_front);
    EXPECT_NE(free_ends_at::none, get<0>(cfg).value & free_ends_at::seq1_back);
    EXPECT_EQ(free_ends_at::none, get<0>(cfg).value & free_ends_at::seq2_back);

    EXPECT_TRUE((std::is_same_v<remove_cvref_t<decltype(cfg)>,
                                detail::configuration<detail::align_config_sequence_ends_deferred>>));
}

TEST(align_config_sequence_ends, invoke_static)
{
    auto cfg = std::invoke(align_cfg::sequence_ends<free_ends_at::seq1>(), detail::configuration<>{});

    EXPECT_NE(free_ends_at::none, get<0>(cfg).value & free_ends_at::seq1_front);
    EXPECT_EQ(free_ends_at::none, get<0>(cfg).value & free_ends_at::seq2_front);
    EXPECT_NE(free_ends_at::none, get<0>(cfg).value & free_ends_at::seq1_back);
    EXPECT_EQ(free_ends_at::none, get<0>(cfg).value & free_ends_at::seq2_back);

    EXPECT_TRUE((std::is_same_v<remove_cvref_t<decltype(cfg)>,
                                detail::configuration<detail::align_config_sequence_ends<free_ends_at::seq1>>>));
}

TEST(align_config_sequence_ends, get_by_enum)
{
    {
        detail::configuration cfg = align_cfg::sequence_ends<>(free_ends_at::seq1_back | free_ends_at::seq2_front);

        EXPECT_EQ(free_ends_at::none, get<0>(cfg).value & free_ends_at::seq1_front);
        EXPECT_NE(free_ends_at::none, get<0>(cfg).value & free_ends_at::seq2_front);
        EXPECT_NE(free_ends_at::none, get<0>(cfg).value & free_ends_at::seq1_back);
        EXPECT_EQ(free_ends_at::none, get<0>(cfg).value & free_ends_at::seq2_back);

        EXPECT_TRUE((std::is_same_v<decltype(get<align_cfg::id::sequence_ends>(cfg)), free_ends_at &>));
    }

    {
        detail::configuration<detail::align_config_sequence_ends_deferred> const c_cfg =
            detail::configuration{align_cfg::sequence_ends<>(free_ends_at::seq1_back | free_ends_at::seq2_front)};

        EXPECT_EQ(free_ends_at::none, get<0>(c_cfg).value & free_ends_at::seq1_front);
        EXPECT_NE(free_ends_at::none, get<0>(c_cfg).value & free_ends_at::seq2_front);
        EXPECT_NE(free_ends_at::none, get<0>(c_cfg).value & free_ends_at::seq1_back);
        EXPECT_EQ(free_ends_at::none, get<0>(c_cfg).value & free_ends_at::seq2_back);

        EXPECT_TRUE((std::is_same_v<decltype(get<align_cfg::id::sequence_ends>(c_cfg)), free_ends_at const &>));
    }

    {
        detail::configuration cfg = align_cfg::sequence_ends<>(free_ends_at::seq1_back | free_ends_at::seq2_front);

        EXPECT_EQ(free_ends_at::none, get<0>(cfg).value & free_ends_at::seq1_front);
        EXPECT_NE(free_ends_at::none, get<0>(cfg).value & free_ends_at::seq2_front);
        EXPECT_NE(free_ends_at::none, get<0>(cfg).value & free_ends_at::seq1_back);
        EXPECT_EQ(free_ends_at::none, get<0>(cfg).value & free_ends_at::seq2_back);

        EXPECT_TRUE((std::is_same_v<decltype(get<align_cfg::id::sequence_ends>(std::move(cfg))), free_ends_at &&>));
    }

    {
        detail::configuration<detail::align_config_sequence_ends_deferred> const c_cfg =
            detail::configuration{align_cfg::sequence_ends<>(free_ends_at::seq1_back | free_ends_at::seq2_front)};

        EXPECT_EQ(free_ends_at::none, get<0>(c_cfg).value & free_ends_at::seq1_front);
        EXPECT_NE(free_ends_at::none, get<0>(c_cfg).value & free_ends_at::seq2_front);
        EXPECT_NE(free_ends_at::none, get<0>(c_cfg).value & free_ends_at::seq1_back);
        EXPECT_EQ(free_ends_at::none, get<0>(c_cfg).value & free_ends_at::seq2_back);

        EXPECT_TRUE((std::is_same_v<decltype(get<align_cfg::id::sequence_ends>(std::move(c_cfg))),
                                    free_ends_at const &&>));
    }
}

TEST(align_config_sequence_ends, free_ends_enum_all_and_none)
{
    seqan3::detail::align_config_sequence_ends<free_ends_at::all> cfg_all;
    EXPECT_NE(free_ends_at::none, cfg_all.value & free_ends_at::seq1_front);
    EXPECT_NE(free_ends_at::none, cfg_all.value & free_ends_at::seq2_front);
    EXPECT_NE(free_ends_at::none, cfg_all.value & free_ends_at::seq1_back);
    EXPECT_NE(free_ends_at::none, cfg_all.value & free_ends_at::seq2_back);

    seqan3::detail::align_config_sequence_ends cfg_none;
    EXPECT_EQ(free_ends_at::none, cfg_none.value & free_ends_at::seq1_front);
    EXPECT_EQ(free_ends_at::none, cfg_none.value & free_ends_at::seq2_front);
    EXPECT_EQ(free_ends_at::none, cfg_none.value & free_ends_at::seq1_back);
    EXPECT_EQ(free_ends_at::none, cfg_none.value & free_ends_at::seq2_back);
}

TEST(align_config_sequence_ends, invoke_deferred)
{
    detail::configuration cfg = align_cfg::sequence_ends<>(free_ends_at::seq1);

    auto call_on_site = [] (auto && new_cfg)
    {
        return std::get<0>(new_cfg).value;
    };

    EXPECT_EQ((std::invoke(std::get<0>(cfg), call_on_site, cfg)), free_ends_at::seq1);
}

// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <functional>
#include <type_traits>

#include <seqan3/alignment/configuration/align_config_score.hpp>
#include <seqan3/alignment/scoring/aminoacid_scoring_scheme.hpp>

using namespace seqan3;

struct bar
{
    int value;
};

TEST(align_config_score, constructor)
{
    EXPECT_TRUE((std::is_default_constructible_v<detail::align_config_score<aminoacid_scoring_scheme<>>>));
}

TEST(align_config_score, on_align_config)
{
    using config_score_type = detail::align_config_score<aminoacid_scoring_scheme<>>;

    EXPECT_TRUE((std::is_same_v<typename detail::on_align_config<align_cfg::id::score>::invoke<config_score_type>,
                 std::true_type>));
    EXPECT_TRUE((std::is_same_v<typename detail::on_align_config<align_cfg::id::score>::invoke<bar>,
                 std::false_type>));
}

TEST(align_config_score, align_config_type_to_id)
{
    using config_score_type = detail::align_config_score<aminoacid_scoring_scheme<>>;

    EXPECT_EQ(detail::align_config_type_to_id<config_score_type>::value, align_cfg::id::score);
    EXPECT_EQ(detail::align_config_type_to_id_v<config_score_type>, align_cfg::id::score);
}

TEST(align_config_score, invoke)
{
    auto cfg = std::invoke(align_cfg::score(aminoacid_scoring_scheme(aminoacid_similarity_matrix::BLOSUM62)),
                           detail::configuration<>{});

    EXPECT_EQ(std::get<0>(cfg).value.score('I'_aa27, 'V'_aa27), 3);
    EXPECT_TRUE((std::is_same_v<remove_cvref_t<decltype(cfg)>,
                                detail::configuration<detail::align_config_score<aminoacid_scoring_scheme<int8_t>>>>));
}

TEST(align_config_score, get_by_enum)
{
    aminoacid_scoring_scheme scheme(aminoacid_similarity_matrix::BLOSUM62);
    detail::configuration cfg = align_cfg::score(scheme);

    EXPECT_EQ(get<align_cfg::id::score>(cfg).score('I'_aa27, 'V'_aa27), 3);
    EXPECT_TRUE((std::is_same_v<decltype(get<align_cfg::id::score>(cfg)),
                                aminoacid_scoring_scheme<int8_t> &>));

    EXPECT_EQ(get<align_cfg::id::score>(std::move(cfg)).score('I'_aa27, 'V'_aa27), 3);
    EXPECT_TRUE((std::is_same_v<decltype(get<align_cfg::id::score>(std::move(cfg))),
                                aminoacid_scoring_scheme<int8_t> &&>));

    detail::configuration<detail::align_config_score<aminoacid_scoring_scheme<>>> const c_cfg =
        detail::configuration{align_cfg::score(scheme)};

    EXPECT_EQ(get<align_cfg::id::score>(c_cfg).score('I'_aa27, 'V'_aa27), 3);
    EXPECT_TRUE((std::is_same_v<decltype(get<align_cfg::id::score>(c_cfg)),
                                aminoacid_scoring_scheme<int8_t> const &>));

    EXPECT_EQ(get<align_cfg::id::score>(std::move(c_cfg)).score('I'_aa27, 'V'_aa27), 3);
    EXPECT_TRUE((std::is_same_v<decltype(get<align_cfg::id::score>(std::move(c_cfg))),
                                aminoacid_scoring_scheme<int8_t> const &&>));
}

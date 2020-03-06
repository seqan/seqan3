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

struct bar
{
    int value;
};

TEST(align_config_score, constructor)
{
    EXPECT_TRUE((std::is_default_constructible_v<seqan3::detail::align_config_score<seqan3::aminoacid_scoring_scheme<>>>));
}

TEST(align_config_score, on_align_config)
{
    using config_score_type = seqan3::detail::align_config_score<seqan3::aminoacid_scoring_scheme<>>;
    using typename_of_align_config = typename seqan3::detail::on_align_config<seqan3::align_cfg::id::score;

    EXPECT_TRUE((std::is_same_v<typename_of_align_config>::invoke<config_score_type>,
                 std::true_type>));
    EXPECT_TRUE((std::is_same_v<typename_of_align_config>::invoke<bar>,
                 std::false_type>));
}

TEST(align_config_score, align_config_type_to_id)
{
    using config_score_type = seqan3::detail::align_config_score<seqan3::aminoacid_scoring_scheme<>>;

    EXPECT_EQ(seqan3::detail::align_config_type_to_id<config_score_type>::value, seqan3::align_cfg::id::score);
    EXPECT_EQ(seqan3::detail::align_config_type_to_id_v<config_score_type>, seqan3::align_cfg::id::score);
}

TEST(align_config_score, invoke)
{
    using blossum_matrix = aminoacid_similarity_matrix::BLOSUM62;
    auto cfg = std::invoke(seqan3::align_cfg::score(seqan3::aminoacid_scoring_scheme(blossum_matrix)),
                           seqan3::detail::configuration<>{});

    EXPECT_EQ(std::get<0>(cfg).value.score('I'_aa27, 'V'_aa27), 3);

    using aminoacid_scoring_scheme_cfg = seqan3::detail::align_config_score<seqan3::aminoacid_scoring_scheme<int8_t>>;
    EXPECT_TRUE((std::is_same_v<remove_cvref_t<decltype(cfg)>,
                                seqan3::detail::configuration<aminoacid_scoring_scheme_cfg>>));
}

TEST(align_config_score, get_by_enum)
{
    seqan3::aminoacid_scoring_scheme scheme(aminoacid_similarity_matrix::BLOSUM62);
    seqan3::detail::configuration cfg = seqan3::align_cfg::score(scheme);

    EXPECT_EQ(std::get<seqan3::align_cfg::id::score>(cfg).score('I'_aa27, 'V'_aa27), 3);
    EXPECT_TRUE((std::is_same_v<decltype(std::get<seqan3::align_cfg::id::score>(cfg)),
                                seqan3::aminoacid_scoring_scheme<int8_t> &>));

    EXPECT_EQ(std::get<seqan3::align_cfg::id::score>(std::move(cfg)).score('I'_aa27, 'V'_aa27), 3);
    EXPECT_TRUE((std::is_same_v<decltype(std::get<seqan3::align_cfg::id::score>(std::move(cfg))),
                                seqan3::aminoacid_scoring_scheme<int8_t> &&>));

    seqan3::detail::configuration<seqan3::detail::align_config_score<seqan3::aminoacid_scoring_scheme<>>> const c_cfg =
        seqan3::detail::configuration{seqan3::align_cfg::score(scheme)};

    EXPECT_EQ(std::get<seqan3::align_cfg::id::score>(c_cfg).score('I'_aa27, 'V'_aa27), 3);
    EXPECT_TRUE((std::is_same_v<decltype(std::get<seqan3::align_cfg::id::score>(c_cfg)),
                                seqan3::aminoacid_scoring_scheme<int8_t> const &>));

    EXPECT_EQ(std::get<seqan3::align_cfg::id::score>(std::move(c_cfg)).score('I'_aa27, 'V'_aa27), 3);
    EXPECT_TRUE((std::is_same_v<decltype(std::get<seqan3::align_cfg::id::score>(std::move(c_cfg))),
                                seqan3::aminoacid_scoring_scheme<int8_t> const &&>));
}

// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ============================================================================

/* \file
 * \brief Provides tests for seqan3::detail::align_config_score.
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 */

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

    EXPECT_EQ(std::get<0>(cfg).value.score(aa27::I, aa27::V), 3);
    EXPECT_TRUE((std::is_same_v<remove_cvref_t<decltype(cfg)>,
                                detail::configuration<detail::align_config_score<aminoacid_scoring_scheme<int8_t>>>>));
}

TEST(align_config_score, get_by_enum)
{
    aminoacid_scoring_scheme scheme(aminoacid_similarity_matrix::BLOSUM62);
    detail::configuration cfg = align_cfg::score(scheme);

    EXPECT_EQ(get<align_cfg::id::score>(cfg).score(aa27::I, aa27::V), 3);
    EXPECT_TRUE((std::is_same_v<decltype(get<align_cfg::id::score>(cfg)),
                                aminoacid_scoring_scheme<int8_t> &>));

    EXPECT_EQ(get<align_cfg::id::score>(std::move(cfg)).score(aa27::I, aa27::V), 3);
    EXPECT_TRUE((std::is_same_v<decltype(get<align_cfg::id::score>(std::move(cfg))),
                                aminoacid_scoring_scheme<int8_t> &&>));

    detail::configuration<detail::align_config_score<aminoacid_scoring_scheme<>>> const c_cfg =
        detail::configuration{align_cfg::score(scheme)};

    EXPECT_EQ(get<align_cfg::id::score>(c_cfg).score(aa27::I, aa27::V), 3);
    EXPECT_TRUE((std::is_same_v<decltype(get<align_cfg::id::score>(c_cfg)),
                                aminoacid_scoring_scheme<int8_t> const &>));

    EXPECT_EQ(get<align_cfg::id::score>(std::move(c_cfg)).score(aa27::I, aa27::V), 3);
    EXPECT_TRUE((std::is_same_v<decltype(get<align_cfg::id::score>(std::move(c_cfg))),
                                aminoacid_scoring_scheme<int8_t> const &&>));
}

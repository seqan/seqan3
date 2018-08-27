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

/*!\file
 * \brief Provides tests for seqan3::align_config_gap.
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 */

#include <gtest/gtest.h>

#include <functional>
#include <type_traits>

#include <seqan3/alignment/configuration/align_config_gap.hpp>
#include <seqan3/alignment/gap/affine.hpp>
#include <seqan3/alignment/gap/linear.hpp>
#include <seqan3/core/detail/reflection.hpp>

using namespace seqan3;

struct bar
{
    int value;
};

TEST(align_config_gap, constructor)
{
    EXPECT_TRUE((std::is_default_constructible_v<detail::align_config_gap<gap_linear<int32_t>>>));
    EXPECT_TRUE((std::is_default_constructible_v<detail::align_config_gap<gap_affine<int32_t>>>));
}

TEST(align_config_gap, on_align_config)
{
    using gap_linear_config_t = detail::align_config_gap<gap_linear<int32_t>>;
    using gap_affine_config_t = detail::align_config_gap<gap_affine<int32_t>>;

    EXPECT_TRUE((std::is_same_v<typename detail::on_align_config<align_cfg::id::gap>::invoke<gap_linear_config_t>,
                 std::true_type>));
    EXPECT_TRUE((std::is_same_v<typename detail::on_align_config<align_cfg::id::gap>::invoke<gap_affine_config_t>,
                 std::true_type>));
    EXPECT_TRUE((std::is_same_v<typename detail::on_align_config<align_cfg::id::gap>::invoke<bar>,
                 std::false_type>));
}

TEST(align_config_gap, align_config_type_to_id)
{
    using gap_linear_config_t = detail::align_config_gap<gap_linear<int32_t>>;
    using gap_affine_config_t = detail::align_config_gap<gap_affine<int32_t>>;

    EXPECT_EQ(detail::align_config_type_to_id<gap_linear_config_t>::value, align_cfg::id::gap);
    EXPECT_EQ(detail::align_config_type_to_id<gap_affine_config_t>::value, align_cfg::id::gap);
    EXPECT_EQ(detail::align_config_type_to_id_v<gap_linear_config_t>, align_cfg::id::gap);
    EXPECT_EQ(detail::align_config_type_to_id_v<gap_affine_config_t>, align_cfg::id::gap);
}

TEST(align_config_gap, invoke_linear)
{
    {
        auto cfg = std::invoke(align_cfg::gap_linear(gap_score{-4}), detail::configuration<>{});

        EXPECT_EQ(std::get<0>(cfg).value.scheme.get_gap_score(), -4);
        EXPECT_TRUE((std::is_same_v<remove_cvref_t<decltype(cfg)>,
                                    detail::configuration<detail::align_config_gap<gap_linear<int32_t>>>>));
    }

    {
        detail::configuration cfg = align_cfg::gap_linear(gap_score{-4.1});

        EXPECT_DOUBLE_EQ(std::get<0>(cfg).value.scheme.get_gap_score(), -4.1);
        EXPECT_TRUE((std::is_same_v<remove_cvref_t<decltype(cfg)>,
                                    detail::configuration<detail::align_config_gap<gap_linear<double>>>>));
    }
}

TEST(align_config_gap, invoke_affine)
{
    {
        auto cfg = std::invoke(align_cfg::gap_affine(gap_score{-2}, gap_open_score{-4}), detail::configuration<>{});

        EXPECT_EQ(get<0>(cfg).value.scheme.get_gap_score(), -2);
        EXPECT_EQ(get<0>(cfg).value.scheme.get_gap_open_score(), -4);
        EXPECT_TRUE((std::is_same_v<remove_cvref_t<decltype(cfg)>,
                                    detail::configuration<detail::align_config_gap<gap_affine<int32_t>>>>));
    }

    {
        detail::configuration cfg = align_cfg::gap_affine(gap_score{-2.2}, gap_open_score{-4.1});

        EXPECT_DOUBLE_EQ(get<0>(cfg).value.scheme.get_gap_score(), -2.2);
        EXPECT_DOUBLE_EQ(get<0>(cfg).value.scheme.get_gap_open_score(), -4.1);
        EXPECT_TRUE((std::is_same_v<remove_cvref_t<decltype(cfg)>,
                                    detail::configuration<detail::align_config_gap<gap_affine<double>>>>));
    }
}

TEST(align_config_gap, get_by_enum_linear)
{
    {
        detail::configuration cfg = align_cfg::gap_linear(gap_score{-4});

        EXPECT_EQ(get<align_cfg::id::gap>(cfg).scheme.get_gap_score(), -4);
        EXPECT_TRUE((std::is_same_v<decltype(get<align_cfg::id::gap>(cfg)),
                                    gap_linear<int> &>));
    }

    { // TODO(rrahn): Investigate why not deduction guide. See #314
        detail::configuration<detail::align_config_gap<gap_linear<int>>> const c_cfg =
            detail::configuration{align_cfg::gap_linear(gap_score{-4})};

        EXPECT_EQ(get<align_cfg::id::gap>(c_cfg).scheme.get_gap_score(), -4);
        EXPECT_TRUE((std::is_same_v<decltype(get<align_cfg::id::gap>(c_cfg)),
                                    gap_linear<int> const &>));
    }

    {
        detail::configuration cfg = align_cfg::gap_linear(gap_score{-4});

        EXPECT_EQ(get<align_cfg::id::gap>(std::move(cfg)).scheme.get_gap_score(), -4);
        EXPECT_TRUE((std::is_same_v<decltype(get<align_cfg::id::gap>(std::move(cfg))),
                                    gap_linear<int> &&>));
    }

    {
        detail::configuration<detail::align_config_gap<gap_linear<int>>> const c_cfg =
            detail::configuration{align_cfg::gap_linear(gap_score{-4})};

        EXPECT_EQ(get<align_cfg::id::gap>(std::move(c_cfg)).scheme.get_gap_score(), -4);
        EXPECT_TRUE((std::is_same_v<decltype(get<align_cfg::id::gap>(std::move(c_cfg))),
                                    gap_linear<int> const &&>));
    }
}

TEST(align_config_gap, get_by_enum_affine)
{
    {
        detail::configuration cfg = align_cfg::gap_affine(gap_score{-2}, gap_open_score{-4});

        EXPECT_EQ(get<align_cfg::id::gap>(cfg).scheme.get_gap_score(), -2);
        EXPECT_EQ(get<align_cfg::id::gap>(cfg).scheme.get_gap_open_score(), -4);
        EXPECT_TRUE((std::is_same_v<decltype(get<align_cfg::id::gap>(cfg)),
                                    gap_affine<int> &>));
    }

    {
        detail::configuration<detail::align_config_gap<gap_affine<int>>> const c_cfg =
            detail::configuration{align_cfg::gap_affine(gap_score{-2}, gap_open_score{-4})};

        EXPECT_EQ(get<align_cfg::id::gap>(c_cfg).scheme.get_gap_score(), -2);
        EXPECT_EQ(get<align_cfg::id::gap>(c_cfg).scheme.get_gap_open_score(), -4);
        EXPECT_TRUE((std::is_same_v<decltype(get<align_cfg::id::gap>(c_cfg)),
                                    gap_affine<int> const &>));
    }

    {
        detail::configuration cfg = align_cfg::gap_affine(gap_score{-2}, gap_open_score{-4});

        EXPECT_EQ(get<align_cfg::id::gap>(std::move(cfg)).scheme.get_gap_score(), -2);
        EXPECT_EQ(get<align_cfg::id::gap>(std::move(cfg)).scheme.get_gap_open_score(), -4);
        EXPECT_TRUE((std::is_same_v<decltype(get<align_cfg::id::gap>(std::move(cfg))),
                                    gap_affine<int> &&>));
    }

    {
        detail::configuration<detail::align_config_gap<gap_affine<int>>> const c_cfg =
            detail::configuration{align_cfg::gap_affine(gap_score{-2}, gap_open_score{-4})};

        EXPECT_EQ(get<align_cfg::id::gap>(std::move(c_cfg)).scheme.get_gap_score(), -2);
        EXPECT_EQ(get<align_cfg::id::gap>(std::move(c_cfg)).scheme.get_gap_open_score(), -4);
        EXPECT_TRUE((std::is_same_v<decltype(get<align_cfg::id::gap>(std::move(c_cfg))),
                                    gap_affine<int> const &&>));
    }
}

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
#include <seqan3/core/detail/reflection.hpp>

using namespace seqan3;

struct bar
{
    int value;
};

TEST(align_config_gap, constructor)
{
    EXPECT_TRUE((std::is_default_constructible_v<detail::align_config_gap<int32_t>>));
}

TEST(align_config_gap, on_align_config)
{
    using gap_config_type = detail::align_config_gap<int32_t>;

    EXPECT_TRUE((std::is_same_v<typename detail::on_align_config<align_cfg::id::gap>::invoke<gap_config_type>,
                 std::true_type>));
    EXPECT_TRUE((std::is_same_v<typename detail::on_align_config<align_cfg::id::gap>::invoke<bar>,
                 std::false_type>));
}

TEST(align_config_gap, align_config_type_to_id)
{
    using gap_config_type = detail::align_config_gap<int32_t>;

    EXPECT_EQ(detail::align_config_type_to_id<gap_config_type>::value, align_cfg::id::gap);
    EXPECT_EQ(detail::align_config_type_to_id_v<gap_config_type>, align_cfg::id::gap);
}

TEST(align_config_gap, invoke_linear)
{
    { // score type: int
        auto cfg = std::invoke(align_cfg::gap(gap_score{-4}), detail::configuration<>{});

        EXPECT_EQ(std::get<0>(cfg).value.get_gap_score(), -4);
        EXPECT_TRUE((std::is_same_v<remove_cvref_t<decltype(cfg)>,
                                    detail::configuration<detail::align_config_gap<gap_scheme<int32_t>>>>));
    }

    { // score type: double
        detail::configuration cfg = align_cfg::gap(gap_score{-4.1});

        EXPECT_DOUBLE_EQ(std::get<0>(cfg).value.get_gap_score(), -4.1);
        EXPECT_TRUE((std::is_same_v<remove_cvref_t<decltype(cfg)>,
                                    detail::configuration<detail::align_config_gap<gap_scheme<double>>>>));
    }

    { // gap_scheme<int> explicit
        gap_scheme<int> scheme{gap_score{-4}};
        detail::configuration cfg = align_cfg::gap(scheme);

        EXPECT_EQ(std::get<0>(cfg).value.get_gap_score(), -4);
        EXPECT_TRUE((std::is_same_v<remove_cvref_t<decltype(cfg)>,
                                    detail::configuration<detail::align_config_gap<gap_scheme<int32_t>>>>));
    }

    { // gap_scheme<double> explicit
        gap_scheme<double> scheme{gap_score{-4.1}};
        detail::configuration cfg = align_cfg::gap(scheme);

        EXPECT_DOUBLE_EQ(std::get<0>(cfg).value.get_gap_score(), -4.1);
        EXPECT_TRUE((std::is_same_v<remove_cvref_t<decltype(cfg)>,
                                    detail::configuration<detail::align_config_gap<gap_scheme<double>>>>));
    }
}

TEST(align_config_gap, invoke_affine)
{
    { // score type: int
        auto cfg = std::invoke(align_cfg::gap(gap_score{-2}, gap_open_score{-4}), detail::configuration<>{});

        EXPECT_EQ(get<0>(cfg).value.get_gap_score(), -2);
        EXPECT_EQ(get<0>(cfg).value.get_gap_open_score(), -4);
        EXPECT_TRUE((std::is_same_v<remove_cvref_t<decltype(cfg)>,
                                    detail::configuration<detail::align_config_gap<gap_scheme<int32_t>>>>));
    }

    { // score type: float
        detail::configuration cfg = align_cfg::gap(gap_score{-2.2f}, gap_open_score{-4.1f});

        EXPECT_FLOAT_EQ(get<0>(cfg).value.get_gap_score(), -2.2f);
        EXPECT_FLOAT_EQ(get<0>(cfg).value.get_gap_open_score(), -4.1f);
        EXPECT_TRUE((std::is_same_v<remove_cvref_t<decltype(cfg)>,
                                    detail::configuration<detail::align_config_gap<gap_scheme<float>>>>));
    }

    { // gap_scheme<int8_t>
        gap_scheme scheme{gap_score{-2}, gap_open_score{-4}}; // int is converted to int8_t
        auto cfg = std::invoke(align_cfg::gap(scheme), detail::configuration<>{});

        EXPECT_EQ(get<0>(cfg).value.get_gap_score(), static_cast<int8_t>(-2));
        EXPECT_EQ(get<0>(cfg).value.get_gap_open_score(), static_cast<int8_t>(-4));
        EXPECT_TRUE((std::is_same_v<remove_cvref_t<decltype(cfg)>,
                                    detail::configuration<detail::align_config_gap<gap_scheme<int8_t>>>>));
    }

    { // gap_scheme<float>
        gap_scheme scheme{gap_score{-2.2}, gap_open_score{-4.1}}; // double is converted to float
        detail::configuration cfg = align_cfg::gap(scheme);

        EXPECT_FLOAT_EQ(get<0>(cfg).value.get_gap_score(), -2.2f);
        EXPECT_FLOAT_EQ(get<0>(cfg).value.get_gap_open_score(), -4.1f);
        EXPECT_TRUE((std::is_same_v<remove_cvref_t<decltype(cfg)>,
                                    detail::configuration<detail::align_config_gap<gap_scheme<float>>>>));
    }
}

TEST(align_config_gap, get_by_enum_linear)
{
    {
        detail::configuration cfg = align_cfg::gap(gap_score{-4});

        EXPECT_EQ(get<align_cfg::id::gap>(cfg).get_gap_score(), -4);
        EXPECT_TRUE((std::is_same_v<decltype(get<align_cfg::id::gap>(cfg)),
                                    gap_scheme<int32_t> &>));
    }

    { // TODO(rrahn): Investigate why not deduction guide. See #314
        detail::configuration<detail::align_config_gap<gap_scheme<int>>> const c_cfg =
            detail::configuration{align_cfg::gap(gap_score{-4})};

        EXPECT_EQ(get<align_cfg::id::gap>(c_cfg).get_gap_score(), -4);
        EXPECT_TRUE((std::is_same_v<decltype(get<align_cfg::id::gap>(c_cfg)),
                                    gap_scheme<int32_t> const &>));
    }

    {
        detail::configuration cfg = align_cfg::gap(gap_score{-4});

        EXPECT_EQ(get<align_cfg::id::gap>(std::move(cfg)).get_gap_score(), -4);
        EXPECT_TRUE((std::is_same_v<decltype(get<align_cfg::id::gap>(std::move(cfg))),
                                    gap_scheme<int32_t> &&>));
    }

    {
        // The type of gap_score is not deduced here: must be int8_t, because of explicit configuration type.
        detail::configuration<detail::align_config_gap<gap_scheme<int8_t>>> const c_cfg =
            detail::configuration{align_cfg::gap(gap_score{ static_cast<int8_t>(-4) })};

        EXPECT_EQ(get<align_cfg::id::gap>(std::move(c_cfg)).get_gap_score(), static_cast<int8_t>(-4));
        EXPECT_TRUE((std::is_same_v<decltype(get<align_cfg::id::gap>(std::move(c_cfg))),
                                    gap_scheme<int8_t> const &&>));
    }

    {
        // The same test with gap_scheme<int8_t>.
        gap_scheme scheme{gap_score{-4}};
        detail::configuration<detail::align_config_gap<gap_scheme<int8_t>>> const c_cfg =
            detail::configuration{align_cfg::gap(scheme)};

        EXPECT_EQ(get<align_cfg::id::gap>(std::move(c_cfg)).get_gap_score(), static_cast<int8_t>(-4));
        EXPECT_TRUE((std::is_same_v<decltype(get<align_cfg::id::gap>(std::move(c_cfg))),
                                    gap_scheme<int8_t> const &&>));
    }
}

TEST(align_config_gap, get_by_enum_affine)
{
    {
        detail::configuration cfg = align_cfg::gap(gap_score{-2}, gap_open_score{-4});

        EXPECT_EQ(get<align_cfg::id::gap>(cfg).get_gap_score(), -2);
        EXPECT_EQ(get<align_cfg::id::gap>(cfg).get_gap_open_score(), -4);
        EXPECT_TRUE((std::is_same_v<decltype(get<align_cfg::id::gap>(cfg)),
                                    gap_scheme<int32_t> &>));
    }

    {
        detail::configuration<detail::align_config_gap<gap_scheme<int>>> const c_cfg =
            detail::configuration{align_cfg::gap(gap_score{-2}, gap_open_score{-4})};

        EXPECT_EQ(get<align_cfg::id::gap>(c_cfg).get_gap_score(), -2);
        EXPECT_EQ(get<align_cfg::id::gap>(c_cfg).get_gap_open_score(), -4);
        EXPECT_TRUE((std::is_same_v<decltype(get<align_cfg::id::gap>(c_cfg)),
                                    gap_scheme<int32_t> const &>));
    }

    {
        detail::configuration cfg = align_cfg::gap(gap_score{-2}, gap_open_score{-4});

        EXPECT_EQ(get<align_cfg::id::gap>(std::move(cfg)).get_gap_score(), -2);
        EXPECT_EQ(get<align_cfg::id::gap>(std::move(cfg)).get_gap_open_score(), -4);
        EXPECT_TRUE((std::is_same_v<decltype(get<align_cfg::id::gap>(std::move(cfg))),
                                    gap_scheme<int32_t> &&>));
    }

    {
        detail::configuration<detail::align_config_gap<gap_scheme<int>>> const c_cfg =
            detail::configuration{align_cfg::gap(gap_score{-2}, gap_open_score{-4})};

        EXPECT_EQ(get<align_cfg::id::gap>(std::move(c_cfg)).get_gap_score(), -2);
        EXPECT_EQ(get<align_cfg::id::gap>(std::move(c_cfg)).get_gap_open_score(), -4);
        EXPECT_TRUE((std::is_same_v<decltype(get<align_cfg::id::gap>(std::move(c_cfg))),
                                    gap_scheme<int32_t> const &&>));
    }

    { // The same test with gap_scheme<float>.
        gap_scheme scheme{gap_score{-2.2}, gap_open_score{-4.1}}; // implicitly converted to float
        detail::configuration<detail::align_config_gap<gap_scheme<float>>> const c_cfg =
            detail::configuration{align_cfg::gap(scheme)};

        EXPECT_EQ(get<align_cfg::id::gap>(std::move(c_cfg)).get_gap_score(), -2.2f);
        EXPECT_EQ(get<align_cfg::id::gap>(std::move(c_cfg)).get_gap_open_score(), -4.1f);
        EXPECT_TRUE((std::is_same_v<decltype(get<align_cfg::id::gap>(std::move(c_cfg))),
                                    gap_scheme<float> const &&>));
    }
}

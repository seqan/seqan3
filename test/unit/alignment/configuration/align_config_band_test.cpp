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
 * \brief Provides tests for alignment band configuration.
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 */

#include <gtest/gtest.h>

#include <functional>
#include <type_traits>

#include <seqan3/alignment/configuration/align_config_band.hpp>
#include <seqan3/alignment/band/static.hpp>

using namespace seqan3;

TEST(align_config_band, constructor)
{
    EXPECT_TRUE((std::is_default_constructible_v<detail::align_config_band<band_static<uint32_t>>>));
}

TEST(align_config_band, on_align_config)
{
    struct bar
    {
        int value;
    };

    using band_config_t = detail::align_config_band<band_static<uint32_t>>;

    EXPECT_TRUE((std::is_same_v<typename detail::on_align_config<align_cfg::id::band>::invoke<band_config_t>,
                 std::true_type>));
    EXPECT_TRUE((std::is_same_v<typename detail::on_align_config<align_cfg::id::band>::invoke<bar>,
                 std::false_type>));
}

TEST(align_config_band, align_config_type_to_id)
{
    using band_config_t = detail::align_config_band<band_static<uint32_t>>;
    EXPECT_EQ(detail::align_config_type_to_id<band_config_t>::value, align_cfg::id::band);
    EXPECT_EQ(detail::align_config_type_to_id_v<band_config_t>, align_cfg::id::band);
}

TEST(align_config_band, invoke_unsigned)
{
    auto cfg = std::invoke(align_cfg::band_static(lower_bound{4u}, upper_bound{5u}), detail::configuration<>{});

    EXPECT_EQ(get<0>(cfg).value.lower_bound, 4u);
    EXPECT_EQ(get<0>(cfg).value.upper_bound, 5u);
    EXPECT_TRUE((std::is_same_v<remove_cvref_t<decltype(cfg)>,
                                detail::configuration<detail::align_config_band<band_static<unsigned>>>>));
}

TEST(align_config_band, invoke_int)
{
    auto cfg = std::invoke(align_cfg::band_static(lower_bound{4}, upper_bound{5}), detail::configuration<>{});

    EXPECT_EQ(get<0>(cfg).value.lower_bound, 4u);
    EXPECT_EQ(get<0>(cfg).value.upper_bound, 5u);
    EXPECT_TRUE((std::is_same_v<remove_cvref_t<decltype(cfg)>,
                                detail::configuration<detail::align_config_band<band_static<unsigned>>>>));
}

TEST(align_config_band, get_by_enum)
{
    {
        detail::configuration cfg = align_cfg::band_static(lower_bound{4u}, upper_bound{5u});

        EXPECT_EQ(get<align_cfg::id::band>(cfg).lower_bound, 4u);
        EXPECT_EQ(get<align_cfg::id::band>(cfg).upper_bound, 5u);
        EXPECT_TRUE((std::is_same_v<decltype(get<align_cfg::id::band>(cfg)),
                                    band_static<unsigned> &>));
    }

    {
        detail::configuration<detail::align_config_band<band_static<uint32_t>>> const c_cfg =
            detail::configuration{align_cfg::band_static(lower_bound{4u}, upper_bound{5u})};

        EXPECT_EQ(get<align_cfg::id::band>(c_cfg).lower_bound, 4u);
        EXPECT_EQ(get<align_cfg::id::band>(c_cfg).upper_bound, 5u);
        EXPECT_TRUE((std::is_same_v<decltype(get<align_cfg::id::band>(c_cfg)),
                                    band_static<uint32_t> const &>));
    }

    {
        detail::configuration cfg = align_cfg::band_static(lower_bound{4}, upper_bound{5});

        EXPECT_EQ(get<align_cfg::id::band>(std::move(cfg)).lower_bound, 4u);
        EXPECT_EQ(get<align_cfg::id::band>(std::move(cfg)).upper_bound, 5u);
        EXPECT_TRUE((std::is_same_v<decltype(get<align_cfg::id::band>(std::move(cfg))),
                                    band_static<unsigned> &&>));
    }

    {
        detail::configuration<detail::align_config_band<band_static<uint32_t>>> const c_cfg =
            detail::configuration{align_cfg::band_static(lower_bound{4u}, upper_bound{5u})};

        EXPECT_EQ(get<align_cfg::id::band>(std::move(c_cfg)).lower_bound, 4u);
        EXPECT_EQ(get<align_cfg::id::band>(std::move(c_cfg)).upper_bound, 5u);
        EXPECT_TRUE((std::is_same_v<decltype(get<align_cfg::id::band>(std::move(c_cfg))),
                                    band_static<uint32_t> const &&>));
    }
}

TEST(align_config_band, negative_values)
{
    band_static bs{lower_bound{-2}, upper_bound{2}};
    EXPECT_EQ(bs.lower_bound, 2);
    EXPECT_EQ(bs.upper_bound, 2);

    auto cfg = std::invoke(align_cfg::band_static(lower_bound{-2}, upper_bound{2}), detail::configuration<>{});
    EXPECT_EQ(get<0>(cfg).value.lower_bound, 2);
    EXPECT_EQ(get<0>(cfg).value.upper_bound, 2);
}

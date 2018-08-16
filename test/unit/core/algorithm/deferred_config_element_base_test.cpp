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

#include <gtest/gtest.h>

#include <meta/meta.hpp>

#include <seqan3/core/algorithm/deferred_config_element_base.hpp>
#include <seqan3/core/algorithm/configuration.hpp>

using namespace seqan3;

template <int _d>
struct bar_static
{
    int value{_d};
};

struct bar : public detail::deferred_config_element_base<bar>
{
    template <typename fn_t, typename configuration_t>
    constexpr auto invoke(fn_t && fn, configuration_t && config) const
    {
        if (std::get<bar>(config).value == 1)
            return fn(config.replace_with(*this, bar_static<1>{}));
        else
            return fn(config.replace_with(*this, bar_static<0>{}));
    }

    int value{1};
};

TEST(deferred_config_element_base, concept)
{
    EXPECT_TRUE(detail::deferred_config_element_concept<bar>);
    EXPECT_FALSE(detail::deferred_config_element_concept<int>);
}

TEST(deferred_config_element_base, standard_construction)
{
    EXPECT_TRUE((std::is_default_constructible_v<bar>));
    EXPECT_TRUE((std::is_copy_constructible_v<bar>));
    EXPECT_TRUE((std::is_move_constructible_v<bar>));
    EXPECT_TRUE((std::is_copy_assignable_v<bar>));
    EXPECT_TRUE((std::is_move_assignable_v<bar>));
}

TEST(deferred_config_element_base, invoke)
{
    {
        detail::configuration<bar> cfg{};
        std::get<0>(cfg).value = 3;

        auto call_on_site = [] (auto && new_cfg)
        {
            EXPECT_TRUE((std::is_same_v<std::remove_reference_t<decltype(new_cfg)>,
                                        detail::configuration<bar_static<0>>>));
            return std::get<0>(new_cfg).value;
        };
        EXPECT_EQ((bar{}(call_on_site, cfg)), 0);
    }

    {
        detail::configuration<bar> cfg{};
        std::get<0>(cfg).value = 1;

        auto call_on_site = [] (auto && new_cfg)
        {
            EXPECT_TRUE((std::is_same_v<std::remove_reference_t<decltype(new_cfg)>,
                                        detail::configuration<bar_static<1>>>));
            return std::get<0>(new_cfg).value;
        };
        EXPECT_EQ((bar{}(call_on_site, cfg)), 1);
    }
}

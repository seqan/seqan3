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

#include <functional>
#include <type_traits>

#include <seqan3/alignment/configuration/align_config_max_error.hpp>

using namespace seqan3;

struct bar
{
    int value;
};

TEST(align_config_max_error, constructor)
{
    EXPECT_TRUE((std::is_default_constructible_v<detail::align_config_max_error>));
}

TEST(align_config_max_error, on_align_config)
{
    using global_config_t = detail::align_config_max_error;
    EXPECT_TRUE((std::is_same_v<typename detail::on_align_config<align_cfg::id::max_error>::invoke<global_config_t>,
                 std::true_type>));
    EXPECT_TRUE((std::is_same_v<typename detail::on_align_config<align_cfg::id::max_error>::invoke<bar>,
                 std::false_type>));
}

TEST(align_config_max_error, align_config_type_to_id)
{
    using global_config_t = detail::align_config_max_error;
    EXPECT_EQ(detail::align_config_type_to_id<global_config_t>::value, align_cfg::id::max_error);
    EXPECT_EQ(detail::align_config_type_to_id_v<global_config_t>, align_cfg::id::max_error);
}

TEST(align_config_max_error, invoke)
{
    detail::configuration cfg = align_cfg::max_error(10);

    EXPECT_TRUE((std::is_same_v<remove_cvref_t<decltype(cfg)>,
                                detail::configuration<detail::align_config_max_error>>));
}

TEST(align_config_max_error, get_by_enum)
{
    detail::configuration cfg = align_cfg::max_error(10);
    auto const c_cfg = detail::configuration{align_cfg::max_error(10)};

    EXPECT_EQ(get<align_cfg::id::max_error>(cfg), 10);
    EXPECT_TRUE((std::is_same_v<decltype(get<align_cfg::id::max_error>(cfg)),
                                uint32_t &>));

    EXPECT_TRUE((std::is_same_v<decltype(get<align_cfg::id::max_error>(c_cfg)),
                                uint32_t const &>));

    EXPECT_TRUE((std::is_same_v<decltype(get<align_cfg::id::max_error>(std::move(cfg))),
                                uint32_t &&>));

    EXPECT_TRUE((std::is_same_v<decltype(get<align_cfg::id::max_error>(std::move(c_cfg))),
                                uint32_t const &&>));
}

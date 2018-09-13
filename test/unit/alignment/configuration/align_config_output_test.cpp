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

#include <type_traits>

#include <seqan3/alignment/configuration/align_config_output.hpp>

using namespace seqan3;

struct bar
{
    int value;
};

TEST(align_config_output, constructor)
{
    EXPECT_TRUE((std::is_default_constructible_v<detail::align_config_output<align_result_key::end>>));
}

TEST(align_config_output, on_align_config)
{
    using output_config_t = detail::align_config_output<align_result_key::trace>;
    EXPECT_TRUE((std::is_same_v<typename detail::on_align_config<align_cfg::id::output>::invoke<output_config_t>,
                 std::true_type>));
    EXPECT_TRUE((std::is_same_v<typename detail::on_align_config<align_cfg::id::output>::invoke<bar>,
                 std::false_type>));
}

TEST(align_config_output, align_config_type_to_id)
{
    using output_config_t = detail::align_config_output<align_result_key::begin>;
    EXPECT_EQ(detail::align_config_type_to_id<output_config_t>::value, align_cfg::id::output);
    EXPECT_EQ(detail::align_config_type_to_id_v<output_config_t>, align_cfg::id::output);
}

TEST(align_config_output, invoke)
{
    detail::configuration cfg = align_cfg::output<align_result_key::score>;

    EXPECT_TRUE((std::is_same_v<remove_cvref_t<decltype(cfg)>,
                                detail::configuration<detail::align_config_output<align_result_key::score>>>));
}

TEST(align_config_output, get_by_enum)
{
    detail::configuration cfg = align_cfg::output<align_result_key::score>;
    auto const c_cfg = detail::configuration{align_cfg::output<align_result_key::score>};

    EXPECT_TRUE((std::is_same_v<decltype(get<align_cfg::id::output>(cfg)),
                                align_result_key &>));

    EXPECT_TRUE((std::is_same_v<decltype(get<align_cfg::id::output>(c_cfg)),
                                align_result_key const &>));

    EXPECT_TRUE((std::is_same_v<decltype(get<align_cfg::id::output>(std::move(cfg))),
                                align_result_key &&>));

    EXPECT_TRUE((std::is_same_v<decltype(get<align_cfg::id::output>(std::move(c_cfg))),
                                align_result_key const &&>));
}

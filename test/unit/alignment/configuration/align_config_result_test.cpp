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

#include <seqan3/alignment/configuration/align_config_result.hpp>
#include <seqan3/core/algorithm/concept.hpp>

using namespace seqan3;

TEST(align_config_max_error, config_element_concept)
{
    EXPECT_TRUE((detail::config_element_concept<align_cfg::result<align_result_key::score>>));
}

TEST(align_config_max_error, construction)
{
    EXPECT_TRUE((std::is_default_constructible_v<align_cfg::result<align_result_key::score>>));
    EXPECT_TRUE((std::is_copy_constructible_v<align_cfg::result<align_result_key::score>>));
    EXPECT_TRUE((std::is_move_constructible_v<align_cfg::result<align_result_key::score>>));
    EXPECT_TRUE((std::is_copy_assignable_v<align_cfg::result<align_result_key::score>>));
    EXPECT_TRUE((std::is_move_assignable_v<align_cfg::result<align_result_key::score>>));
}

TEST(align_config_max_error, configuration)
{
    {
        align_cfg::result<align_result_key::score> elem{};
        configuration cfg{elem};
        EXPECT_EQ((std::is_same_v<std::remove_reference_t<
                                    decltype(get<align_cfg::result<align_result_key::score>>(cfg).value)>,
                                  align_result_key>), true);

        EXPECT_EQ((get<align_cfg::result<align_result_key::score>>(cfg).value), align_result_key::score);
    }

    {
        configuration cfg{align_cfg::result<align_result_key::score>{}};
        EXPECT_EQ((std::is_same_v<std::remove_reference_t<
                                    decltype(get<align_cfg::result<align_result_key::score>>(cfg).value)>,
                                  align_result_key>), true);
        EXPECT_EQ((get<align_cfg::result<align_result_key::score>>(cfg).value), align_result_key::score);
    }
}

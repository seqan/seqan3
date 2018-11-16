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

#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include <array>

#include <seqan3/core/algorithm/configuration_utility.hpp>
#include <seqan3/core/algorithm/pipeable_config_element.hpp>

enum class test_algo_id : uint8_t
{
    bar_id = 0,
    bax_id = 1,
    foo_id = 2,
    foobar_id = 3,
    SIZE = 4
};

namespace seqan3::detail
{

template <>
inline constexpr std::array<std::array<bool, static_cast<uint8_t>(::test_algo_id::SIZE)>,
                            static_cast<uint8_t>(::test_algo_id::SIZE)> compatibility_table<::test_algo_id>
{{
    {0, 1, 1, 1},
    {1, 0, 1, 0},
    {1, 1, 0, 1},
    {1, 0, 1, 0}
}};

}  // namespace seqan3::detail

struct bar : public seqan3::pipeable_config_element
{
    static constexpr test_algo_id id{test_algo_id::bar_id};
    int value{1};
};

struct bax : public seqan3::pipeable_config_element
{
    static constexpr test_algo_id id{test_algo_id::bax_id};
    float value{2.2};
};

struct foo : public seqan3::pipeable_config_element
{
    static constexpr test_algo_id id{test_algo_id::foo_id};
    std::string value{"test"};
};

template <typename t = std::vector<int>>
struct foobar : public seqan3::pipeable_config_element
{
    static constexpr test_algo_id id{test_algo_id::foobar_id};
    t value{0, 1, 2, 3};
};

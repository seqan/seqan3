// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2017, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2017, Knut Reinert & MPI Molekulare Genetik
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

#include <sstream>
#include <vector>

#include <gtest/gtest.h>

#include <range/v3/view/zip.hpp>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/range/detail/any_sized_view.hpp>
#include <seqan3/range/view/to_char.hpp>

using namespace seqan3;
using namespace seqan3::literal;

template <typename T>
class any_sized_view : public ::testing::Test
{};

using parameter_types = ::testing::Types<std::integral_constant<ranges::category, ranges::category::input>,
                                         std::integral_constant<ranges::category, ranges::category::forward>,
                                         std::integral_constant<ranges::category, ranges::category::bidirectional>,
                                         std::integral_constant<ranges::category, ranges::category::random_access>>;

TYPED_TEST_CASE(any_sized_view, parameter_types);

TYPED_TEST(any_sized_view, basic)
{
    dna5_vector vec{"ACTTTGATA"_dna5};
    std::string cmp{"ACTTTGATA"};

    seqan3::detail::any_sized_view<char, TypeParam::value> v = vec | view::to_char;

    EXPECT_EQ(ranges::size(v), ranges::size(vec));

    if constexpr(TypeParam::value == ranges::category::random_access)
    {
        for (size_t i = 0; i < ranges::size(v); ++i)
            EXPECT_EQ(v[i], cmp[i]);
    } else
    {
        for (auto [l, r] : ranges::view::zip(cmp, v))
            EXPECT_EQ(l, r);
    }
}

TYPED_TEST(any_sized_view, concepts)
{
    switch (TypeParam::value)
    {
        case ranges::category::random_access:
            EXPECT_TRUE((random_access_range_concept<seqan3::detail::any_sized_view<char &, TypeParam::value>>));
            [[fallthrough]];
        case ranges::category::bidirectional:
            EXPECT_TRUE((bidirectional_range_concept<seqan3::detail::any_sized_view<char &, TypeParam::value>>));
            [[fallthrough]];
        case ranges::category::forward:
            EXPECT_TRUE((forward_range_concept<seqan3::detail::any_sized_view<char &, TypeParam::value>>));
            [[fallthrough]];
        default:
            EXPECT_TRUE((input_range_concept<seqan3::detail::any_sized_view<char &, TypeParam::value>>));
            EXPECT_TRUE((sized_range_concept<seqan3::detail::any_sized_view<char &, TypeParam::value>>));
    }
}

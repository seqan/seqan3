// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
//
// Copyright (c) 2006-2017, Knut Reinert, FU Berlin
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
// ==========================================================================

#include <iostream>

#include <gtest/gtest.h>

#include <range/v3/view/take_exactly.hpp>

#include <seqan3/alphabet/quality/all.hpp>
#include <seqan3/range/view/to_char.hpp>
#include <seqan3/range/view/join.hpp>

using namespace seqan3;
using namespace seqan3::literal;

template <typename T>
class view_join : public ::testing::Test
{};

using view_join_ra_types = ::testing::Types<
    std::integral_constant<seqan3::view_join_flags, seqan3::view_join_flags::DEFAULT>,
    std::integral_constant<seqan3::view_join_flags, seqan3::view_join_flags::SPARSE>,
    std::integral_constant<seqan3::view_join_flags, seqan3::view_join_flags::LAZY>,
    std::integral_constant<seqan3::view_join_flags, seqan3::view_join_flags::SPARSE | seqan3::view_join_flags::LAZY>>;

TYPED_TEST_CASE(view_join, view_join_ra_types);

TYPED_TEST(view_join, view_join_ra__basic)
{
    constexpr bool is_lazy = (TypeParam::value & seqan3::view_join_flags::LAZY) == seqan3::view_join_flags::LAZY;

    std::vector<dna5_vector> vec{"AAAAA"_dna5, "CCCC"_dna5, "GGG"_dna5, "TT"_dna5};
    detail::view_join_ra<std::vector<dna5_vector>, TypeParam::value> v{vec};

    // size()
    if constexpr (!is_lazy)
    {
        EXPECT_EQ(ranges::size(v), 14u);
        EXPECT_FALSE(v.empty());
    }

    // random access
    EXPECT_EQ(v[0], dna5::A);
    EXPECT_EQ(v[5], dna5::C);
    EXPECT_EQ(v[13], dna5::T);

    // begin and end
    EXPECT_EQ(*v.begin(), dna5::A);
    EXPECT_NE(v.begin(), v.end());
    EXPECT_EQ(v.begin() + 14, v.end());

    // front and back
    EXPECT_EQ(v.front(), dna5::A);
    if constexpr (!is_lazy)
    {
        EXPECT_EQ(v.back(), dna5::T);
    }

    // implicit conversion back to container
    dna5_vector cont = v;
    EXPECT_EQ(cont, "AAAAACCCCGGGTT"_dna5);

    // pipable
    EXPECT_EQ(std::string(v | view::to_char), "AAAAACCCCGGGTT");
}

TYPED_TEST(view_join, view_join_ra__concepts)
{
    constexpr bool is_lazy = (TypeParam::value & seqan3::view_join_flags::LAZY) == seqan3::view_join_flags::LAZY;
    using t = detail::view_join_ra<std::vector<dna5_vector>, TypeParam::value>;
    EXPECT_TRUE(input_range_concept<t>);
    EXPECT_TRUE(forward_range_concept<t>);
    EXPECT_TRUE(random_access_range_concept<t>);
    EXPECT_TRUE(view_concept<t>);
    EXPECT_EQ(sized_range_concept<t>, !is_lazy);
}

TYPED_TEST(view_join, join_fn__input_is_ra_range)
{
    std::vector<dna5_vector> vec{"AAAAA"_dna5, "CCCC"_dna5, "GGG"_dna5, "TT"_dna5};

    // pipe notation
    auto v = vec | view::join<TypeParam::value>;
    EXPECT_EQ(std::string(v | view::to_char), "AAAAACCCCGGGTT");

    // function syntax
    auto v2 = view::join<TypeParam::value>(vec);
    EXPECT_EQ(std::string(v2 | view::to_char), "AAAAACCCCGGGTT");

    // combinability
    std::string v3 = vec | view::join<TypeParam::value> | ranges::view::take_exactly(5) | view::to_char;
    EXPECT_EQ("AAAAA", v3);
}

TYPED_TEST(view_join, join_fn__input_is_concatenated_sequences)
{
    concatenated_sequences<dna5_vector> vec{"AAAAA"_dna5, "CCCC"_dna5, "GGG"_dna5, "TT"_dna5};

    // pipe notation
    auto v = vec | view::join<TypeParam::value>;
    EXPECT_EQ(std::string(v | view::to_char), "AAAAACCCCGGGTT");

    // function syntax
    auto v2 = view::join<TypeParam::value>(vec);
    EXPECT_EQ(std::string(v2 | view::to_char), "AAAAACCCCGGGTT");

    // combinability
    std::string v3 = vec | view::join<TypeParam::value> | ranges::view::take_exactly(5) | view::to_char;
    EXPECT_EQ("AAAAA", v3);

    // different type
    EXPECT_TRUE((std::is_same_v<decltype(v2), decltype(vec.concat())>));
}

TYPED_TEST(view_join, join_fn__input_is_input_range)
{
    std::forward_list<dna5_vector> vec{"AAAAA"_dna5, "CCCC"_dna5, "GGG"_dna5, "TT"_dna5};

    // pipe notation
    auto v = vec | view::join<TypeParam::value>;
    EXPECT_EQ(std::string(v | view::to_char), "AAAAACCCCGGGTT");

    // function syntax
    auto v2 = view::join<TypeParam::value>(vec);
    EXPECT_EQ(std::string(v2 | view::to_char), "AAAAACCCCGGGTT");

    // combinability
    std::string v3 = vec | view::join<TypeParam::value> | ranges::view::take_exactly(5) | view::to_char;
    EXPECT_EQ("AAAAA", v3);

    // different type
    EXPECT_TRUE((std::is_same_v<decltype(v2), decltype(ranges::view::join(vec))>));
}

//TODO tests for ranges::view::join forward

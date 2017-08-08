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

#include <gtest/gtest.h>
#include <sstream>

#include <range/v3/range_traits.hpp>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/range/container/all.hpp>
#include <seqan3/range/view/to_char.hpp>

using namespace seqan3;
using namespace seqan3::literal;

template <typename T>
class container_of_container : public ::testing::Test
{};

using container_of_container_types = ::testing::Types<std::vector<std::vector<dna4>>,
                                                      concatenated_sequences<std::vector<dna4>>>;

TYPED_TEST_CASE(container_of_container, container_of_container_types);

TYPED_TEST(container_of_container, concepts)
{
    EXPECT_TRUE(container_concept<TypeParam>);
    EXPECT_TRUE(container_concept<ranges::range_value_type_t<TypeParam>>);
}

TYPED_TEST(container_of_container, construction)
{
    TypeParam t1;
    TypeParam t2{};
//     EXPECT_EQ(t1, t2); // doesnt work because of conflict between ranges and gtest
    EXPECT_TRUE(t1 == t2);

    // initializer list
    TypeParam t3{"ACGT"_dna4, "ACGT"_dna4, "GAGGA"_dna4};
    TypeParam t4 = {"ACGT"_dna4, "ACGT"_dna4, "GAGGA"_dna4};
    EXPECT_TRUE(t3 == t4);

    // n * value
    TypeParam t5{2, "ACGT"_dna4};
    // from another TypeParam's sub-range
    TypeParam t6{t3.begin(), t3.begin() + 2};
    EXPECT_TRUE(t5 == t6);

    std::vector<std::vector<dna4>> other_vector{"ACGT"_dna4, "ACGT"_dna4, "GAGGA"_dna4};
    // direct from another container-of-container
    TypeParam t7{other_vector};
    // from another container-of-container's sub-range
    TypeParam t8{other_vector.cbegin(), other_vector.cend()};
    EXPECT_TRUE(t3 == t7);
    EXPECT_TRUE(t7 == t8);
}

TYPED_TEST(container_of_container, assign)
{
    TypeParam t1{"ACGT"_dna4, "ACGT"_dna4, "GAGGA"_dna4};
    TypeParam t2{"ACGT"_dna4, "ACGT"_dna4};
    std::vector<std::vector<dna4>> other_vector{"ACGT"_dna4, "ACGT"_dna4, "GAGGA"_dna4};

    // n * value
    TypeParam t3;
    t3.assign(2, "ACGT"_dna4);
    EXPECT_TRUE(t3 == t2);

    // from another container-of-container's sub-range
    TypeParam t4;
    t4.assign(other_vector.cbegin(), other_vector.cend());
    EXPECT_TRUE(t4 == t1);

    // initializer list
    TypeParam t5, t6;
    t5.assign({"ACGT"_dna4, "ACGT"_dna4, "GAGGA"_dna4});
    t6 = {"ACGT"_dna4, "ACGT"_dna4, "GAGGA"_dna4};
    EXPECT_TRUE(t5 == t1);
    EXPECT_TRUE(t6 == t1);

    // direct from another container-of-container
    if constexpr (std::is_same_v<TypeParam, concatenated_sequences<std::vector<dna4>>>)
    {
        TypeParam t7, t8;
        t7.assign(other_vector);
        t8 = other_vector;
        EXPECT_TRUE(t7 == t1);
        EXPECT_TRUE(t8 == t1);
    }
}

TYPED_TEST(container_of_container, iterators)
{
    TypeParam t1{"ACGT"_dna4, "ACGT"_dna4, "GAGGA"_dna4};
    TypeParam const t2{"ACGT"_dna4, "ACGT"_dna4, "GAGGA"_dna4};

    // begin
    EXPECT_TRUE(dna4_vector(*t1.begin())  == "ACGT"_dna4);
    EXPECT_TRUE(dna4_vector(*t1.cbegin()) == "ACGT"_dna4);
    EXPECT_TRUE(dna4_vector(*t2.begin())  == "ACGT"_dna4);
    EXPECT_TRUE(dna4_vector(*t2.cbegin()) == "ACGT"_dna4);

    // end and arithmetic
    EXPECT_TRUE(dna4_vector(*(t1.end()  - 1)) == "GAGGA"_dna4);
    EXPECT_TRUE(dna4_vector(*(t1.cend() - 1)) == "GAGGA"_dna4);
    EXPECT_TRUE(dna4_vector(*(t2.end()  - 1)) == "GAGGA"_dna4);
    EXPECT_TRUE(dna4_vector(*(t2.cend() - 1)) == "GAGGA"_dna4);

    // convertibility between const and non-const
    EXPECT_TRUE(t1.cend() == t1.end());

    // writability
    (*t1.begin())[0] = dna4::T;
    EXPECT_TRUE(dna4_vector(*t1.begin())  == "TCGT"_dna4);
}

TYPED_TEST(container_of_container, element_access)
{
    TypeParam t1{"ACGT"_dna4, "ACGT"_dna4, "GAGGA"_dna4};
    TypeParam const t2{"ACGT"_dna4, "ACGT"_dna4, "GAGGA"_dna4};

    // at
    EXPECT_TRUE(dna4_vector(t1.at(0)) == "ACGT"_dna4);
    EXPECT_TRUE(dna4_vector(t2.at(0)) == "ACGT"_dna4);
    //TODO once we have throwing assert, check at's ability to throw

    // []
    EXPECT_TRUE(dna4_vector(t1[0])    == "ACGT"_dna4);
    EXPECT_TRUE(dna4_vector(t2[0])    == "ACGT"_dna4);

    // front
    EXPECT_TRUE(dna4_vector(t1.front()) == "ACGT"_dna4);
    EXPECT_TRUE(dna4_vector(t2.front()) == "ACGT"_dna4);

    // back
    EXPECT_TRUE(dna4_vector(t1.back()) == "GAGGA"_dna4);
    EXPECT_TRUE(dna4_vector(t2.back()) == "GAGGA"_dna4);

    if constexpr (std::is_same_v<TypeParam, concatenated_sequences<std::vector<dna4>>>)
    {
        using size_type = typename TypeParam::size_type;
        // concat
        EXPECT_TRUE(dna4_vector(t1.concat()) == "ACGTACGTGAGGA"_dna4);
        EXPECT_TRUE(dna4_vector(t2.concat()) == "ACGTACGTGAGGA"_dna4);

        // data
        EXPECT_TRUE(std::get<0>(t1.data()) == "ACGTACGTGAGGA"_dna4);
        EXPECT_TRUE(std::get<0>(t2.data()) == "ACGTACGTGAGGA"_dna4);
        EXPECT_TRUE((std::get<1>(t1.data()) == std::vector<size_type>{0, 4, 8, 13}));
        EXPECT_TRUE((std::get<1>(t2.data()) == std::vector<size_type>{0, 4, 8, 13}));
    }

}

TYPED_TEST(container_of_container, capacity)
{
    TypeParam t0{};
    TypeParam t1{"ACGT"_dna4, "ACGT"_dna4, "GAGGA"_dna4};
    TypeParam const t2{"ACGT"_dna4, "ACGT"_dna4, "GAGGA"_dna4};

    // empty
    EXPECT_TRUE(t0.empty());
    EXPECT_FALSE(t1.empty());
    EXPECT_FALSE(t2.empty());

    // size
    EXPECT_EQ(t0.size(), 0u);
    EXPECT_EQ(t1.size(), 3u);
    EXPECT_EQ(t2.size(), 3u);

    // max_size
    EXPECT_GT(t0.max_size(), 1'000'000'000'000u);
    EXPECT_GT(t1.max_size(), 1'000'000'000'000u);
    EXPECT_GT(t2.max_size(), 1'000'000'000'000u);

    // capacity
    EXPECT_GE(t0.capacity(), t0.size());
    EXPECT_GE(t1.capacity(), t1.size());
    EXPECT_GE(t2.capacity(), t2.size());

    // reserve
    EXPECT_LT(t0.capacity(), 1000u);
    t0.reserve(1000);
    EXPECT_GE(t0.capacity(), 1000u);

    // shrink_to_fit
    t1.reserve(1000);
    EXPECT_GT(t1.capacity(), t1.size()*2);
    t1.shrink_to_fit();
    EXPECT_LE(t1.capacity(), t1.size()*2);

    if constexpr (std::is_same_v<TypeParam, concatenated_sequences<std::vector<dna4>>>)
    {
        // size
        EXPECT_EQ(t0.concat_size(), 0u);
        EXPECT_EQ(t1.concat_size(), 13u);
        EXPECT_EQ(t2.concat_size(), 13u);

        // capacity
        EXPECT_GE(t0.concat_capacity(), t0.concat_size());
        EXPECT_GE(t1.concat_capacity(), t1.concat_size());
        EXPECT_GE(t2.concat_capacity(), t2.concat_size());

        // reserve
        EXPECT_LT(t0.concat_capacity(), 1000u);
        t0.concat_reserve(1000);
        EXPECT_GE(t0.concat_capacity(), 1000u);
    }
}

TYPED_TEST(container_of_container, clear)
{
    TypeParam t0{};
    TypeParam t1{"ACGT"_dna4, "ACGT"_dna4, "GAGGA"_dna4};

    t1.clear();
    EXPECT_TRUE(t0 == t1);
}

TYPED_TEST(container_of_container, insert)
{
    TypeParam t0{};
    TypeParam t1{"ACGT"_dna4, "ACGT"_dna4, "GAGGA"_dna4};

    // position, value
    t0.insert(t0.cend(), "ACGT"_dna4);
    t0.insert(t0.cend(), "GAGGA"_dna4);
    t0.insert(t0.cbegin() + 1, "ACGT"_dna4);
    EXPECT_TRUE(t0 == t1);

    // position, n times values
    t0.clear();
    t1 = {"GAGGA"_dna4, "ACGT"_dna4, "ACGT"_dna4, "GAGGA"_dna4};
    t0.insert(t0.cend(), 2, "ACGT"_dna4);
    t0.insert(t0.cend(), 1, "GAGGA"_dna4);
    t0.insert(t0.cbegin(), 1, "GAGGA"_dna4);
    EXPECT_TRUE(t0 == t1);

    // iterator pair
    t0.clear();
    t1 = {"GAGGA"_dna4, "ACGT"_dna4, "ACGT"_dna4, "GAGGA"_dna4};
    t0.insert(t0.cend(), t1.begin() + 1, t1.begin() + 3);

    t0.insert(t0.cend(),   t1.cend() - 1, t1.cend());
    t0.insert(t0.cbegin(), t1.cend() - 1, t1.cend());
    EXPECT_TRUE(t0 == t1);

    // initializer list
    t0.clear();
    t1 = {"ACGT"_dna4, "ACGT"_dna4, "GAGGA"_dna4};
    t0.insert(t0.cend(), {"ACGT"_dna4, "GAGGA"_dna4});
    t0.insert(t0.cbegin() + 1, "ACGT"_dna4);
    EXPECT_TRUE(t0 == t1);
}

TYPED_TEST(container_of_container, erase)
{
    TypeParam t1{"ACGT"_dna4, "ACGT"_dna4, "GAGGA"_dna4};

    // one element
    t1.erase(t1.begin());
    EXPECT_TRUE((t1 == TypeParam{"ACGT"_dna4, "GAGGA"_dna4}));

    // range
    t1 = {"GAGGA"_dna4, "ACGT"_dna4, "ACGT"_dna4, "GAGGA"_dna4};
    t1.erase(t1.begin() + 1, t1.begin() + 3);
    EXPECT_TRUE((t1 == TypeParam{"GAGGA"_dna4, "GAGGA"_dna4}));
}

TYPED_TEST(container_of_container, push_pop)
{
    TypeParam t0{};

    // push_back
    t0.push_back("ACGT"_dna4);
    EXPECT_TRUE((t0 == TypeParam{"ACGT"_dna4}));
    t0.push_back("GAGGA"_dna4);
    EXPECT_TRUE((t0 == TypeParam{"ACGT"_dna4, "GAGGA"_dna4}));

    // pop_back
    t0.pop_back();
    EXPECT_TRUE(t0 == TypeParam{"ACGT"_dna4});
    t0.pop_back();
    EXPECT_TRUE(t0 == TypeParam{});
}

TYPED_TEST(container_of_container, resize)
{
    TypeParam t0{};

    // enlarge without values
    t0.resize(3);
    EXPECT_TRUE((t0 == TypeParam{{}, {}, {}}));

    // enlarge with value
    t0.resize(5, "ACGT"_dna4);
    EXPECT_TRUE((t0 == TypeParam{{}, {}, {}, "ACGT"_dna4, "ACGT"_dna4}));

    // shrink with value
    t0.resize(4, "ACGT"_dna4);
    EXPECT_TRUE((t0 == TypeParam{{}, {}, {}, "ACGT"_dna4}));

    // shrink without value
    t0.resize(2);
    EXPECT_TRUE((t0 == TypeParam{{}, {}}));
}

TYPED_TEST(container_of_container, swap)
{
    TypeParam t0{};
    TypeParam t1{"ACGT"_dna4, "ACGT"_dna4, "GAGGA"_dna4};

    t0.swap(t1);
    EXPECT_TRUE((t0 == TypeParam{"ACGT"_dna4, "ACGT"_dna4, "GAGGA"_dna4}));
    EXPECT_TRUE((t1 == TypeParam{}));
}

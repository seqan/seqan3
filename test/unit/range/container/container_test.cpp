// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
//
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
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
// ==========================================================================

#include <gtest/gtest.h>

#include <range/v3/algorithm/equal.hpp>
#include <range/v3/range_traits.hpp>
#include <range/v3/view/take.hpp>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/io/stream/debug_stream.hpp>
#include <seqan3/range/container/all.hpp>
#include <seqan3/range/view/convert.hpp>

#if SEQAN3_WITH_CEREAL
#include <seqan3/test/tmp_filename.hpp>

#include <fstream>

#include <cereal/archives/binary.hpp>
#include <cereal/archives/portable_binary.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/xml.hpp>
#include <cereal/types/vector.hpp>
#endif // SEQAN3_WITH_CEREAL

using namespace seqan3;
using namespace seqan3::literal;

template <typename T>
class container : public ::testing::Test
{};

using container_types = ::testing::Types<std::vector<dna4>,
                                         bitcompressed_vector<dna4>>;

TYPED_TEST_CASE(container, container_types);

TYPED_TEST(container, concepts)
{
    EXPECT_TRUE(reservable_container_concept<TypeParam>);
}

TYPED_TEST(container, construction)
{
    TypeParam t1;
    TypeParam t2{};
    EXPECT_EQ(t1, t2);

    // initializer list
    TypeParam t3{dna4::A, dna4::C, dna4::C, dna4::G, dna4::T};
    TypeParam t4 = {dna4::A, dna4::C, dna4::C, dna4::G, dna4::T};
    EXPECT_EQ(t3, t4);

    // n * value
    TypeParam t5{2, dna4::C};

    // from another TypeParam's sub-range
    TypeParam t6{t3.begin() + 1, t3.begin() + 3};
    EXPECT_EQ(t5, t6);

    // direct from another container
    TypeParam t7{"ACCGT"_dna4};
    EXPECT_EQ(t3, t7);
}

TYPED_TEST(container, assign)
{
    TypeParam t0{dna4::C, dna4::C};
    TypeParam t1{dna4::A, dna4::C, dna4::C, dna4::G, dna4::T};

    // n * value
    TypeParam t3;
    t3.assign(2, dna4::C);
    EXPECT_EQ(t3, t0);

    // from another container's sub-range
    TypeParam t4;
    t4.assign(t1.cbegin(), t1.cend());
    EXPECT_EQ(t4, t1);

    // initializer list
    TypeParam t5, t6;
    t5.assign({dna4::A, dna4::C, dna4::C, dna4::G, dna4::T});
    t6 = {dna4::A, dna4::C, dna4::C, dna4::G, dna4::T};
    EXPECT_EQ(t5, t1);
    EXPECT_EQ(t6, t1);

    // direct from another container
    if constexpr (!std::is_same_v<TypeParam, std::vector<dna4>>)
    {
        TypeParam t7;
        t7.assign("ACCGT"_dna4);
        EXPECT_EQ(t7, t1);
    }
}

TYPED_TEST(container, iterators)
{
    TypeParam t1{dna4::A, dna4::C, dna4::C, dna4::G, dna4::T};
    TypeParam const t2{dna4::A, dna4::C, dna4::C, dna4::G, dna4::T};

    // begin
    EXPECT_EQ(*t1.begin(),  dna4::A);
    EXPECT_EQ(*t1.cbegin(), dna4::A);
    EXPECT_EQ(*t2.begin(),  dna4::A);
    EXPECT_EQ(*t2.cbegin(), dna4::A);

    // end and arithmetic
    EXPECT_EQ(*(t1.end()  - 1), dna4::T);
    EXPECT_EQ(*(t1.cend() - 1), dna4::T);
    EXPECT_EQ(*(t2.end()  - 1), dna4::T);
    EXPECT_EQ(*(t2.cend() - 1), dna4::T);

    // convertibility between const and non-const
    EXPECT_TRUE(t1.cend() == t1.end());

    // mutability
    *t1.begin() = dna4::T;
    EXPECT_TRUE((std::ranges::equal(t1, "TCCGT"_dna4)));
}

TYPED_TEST(container, element_access)
{
    TypeParam t1{dna4::A, dna4::C, dna4::C, dna4::G, dna4::T};
    TypeParam const t2{dna4::A, dna4::C, dna4::C, dna4::G, dna4::T};

    // at
    EXPECT_EQ(t1.at(0), dna4::A);
    EXPECT_EQ(t2.at(0), dna4::A);

    //TODO check at's ability to throw

    // []
    EXPECT_EQ(t1[0], dna4::A);
    EXPECT_EQ(t2[0], dna4::A);

    // front
    EXPECT_EQ(t1.front(), dna4::A);
    EXPECT_EQ(t2.front(), dna4::A);

    // back
    EXPECT_EQ(t1.back(), dna4::T);
    EXPECT_EQ(t2.back(), dna4::T);

    // mutability
    t1[0] = dna4::T;
    EXPECT_TRUE((std::ranges::equal(t1, "TCCGT"_dna4)));

    t1.front() = dna4::C;
    EXPECT_TRUE((std::ranges::equal(t1, "CCCGT"_dna4)));

    t1.back() = dna4::G;
    EXPECT_TRUE((std::ranges::equal(t1, "CCCGG"_dna4)));
}

TYPED_TEST(container, capacity)
{
    TypeParam t0{};
    TypeParam t1{dna4::A, dna4::C, dna4::C, dna4::G, dna4::T};
    TypeParam const t2{dna4::A, dna4::C, dna4::C, dna4::G, dna4::T};

    // empty
    EXPECT_TRUE(t0.empty());
    EXPECT_FALSE(t1.empty());
    EXPECT_FALSE(t2.empty());

    // size
    EXPECT_EQ(t0.size(), 0u);
    EXPECT_EQ(t1.size(), 5u);
    EXPECT_EQ(t2.size(), 5u);

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
    EXPECT_LE(t1.capacity(), std::max<size_t>(t1.size()*2, 32ul));
}

TYPED_TEST(container, clear)
{
    TypeParam t0{};
    TypeParam t1{dna4::A, dna4::C, dna4::C, dna4::G, dna4::T};

    t1.clear();
    EXPECT_EQ(t0, t1);
}

TYPED_TEST(container, insert)
{
    TypeParam t0{};
    TypeParam t1{dna4::A, dna4::C, dna4::C, dna4::G, dna4::T};

    // position, value
    t0.insert(t0.cend(), dna4::A);
    t0.insert(t0.cend(), dna4::C);
    t0.insert(t0.cend(), dna4::G);
    t0.insert(t0.cend(), dna4::T);
    t0.insert(t0.cbegin() + 1, dna4::C);
    EXPECT_EQ(t0, t1);

    // position, n times values
    t0.clear();
    t0.insert(t0.cend(), 2, dna4::C);
    t0.insert(t0.cend(), 1, dna4::G);
    t0.insert(t0.cend(), 1, dna4::T);
    t0.insert(t0.cbegin(), 1, dna4::A);
    EXPECT_EQ(t0, t1);

    // iterator pair
    t0.clear();
    t0.insert(t0.cend(), t1.begin() + 1, t1.begin() + 3);

    t0.insert(t0.cend(),   t1.cend() - 2, t1.cend());
    t0.insert(t0.cbegin(), t1.cbegin(), t1.cbegin() + 1);
    EXPECT_EQ(t0, t1);

    // initializer list
    t0.clear();
    t0.insert(t0.cend(), {dna4::A, dna4::C, dna4::G, dna4::T});
    t0.insert(t0.cbegin() + 1, dna4::C);
    EXPECT_EQ(t0, t1);
}

TYPED_TEST(container, erase)
{
    TypeParam t1{dna4::A, dna4::C, dna4::C, dna4::G, dna4::T};

    // one element
    t1.erase(t1.begin());
    EXPECT_EQ(t1, (TypeParam{dna4::C, dna4::C, dna4::G, dna4::T}));

    // range
    t1.erase(t1.begin() + 1, t1.begin() + 3);
    EXPECT_EQ(t1, (TypeParam{dna4::C, dna4::T}));
}

TYPED_TEST(container, push_pop)
{
    TypeParam t0{};

    // push_back
    t0.push_back(dna4::A);
    EXPECT_EQ(t0,  (TypeParam{dna4::A}));
    t0.push_back(dna4::C);
    EXPECT_EQ(t0, (TypeParam{dna4::A, dna4::C}));

    // pop_back
    t0.pop_back();
    EXPECT_EQ(t0, (TypeParam{dna4::A}));
    t0.pop_back();
    EXPECT_EQ(t0, (TypeParam{}));
}

TYPED_TEST(container, resize)
{
    TypeParam t0{};

    // enlarge without values
    t0.resize(3);
    EXPECT_EQ(t0, (TypeParam{dna4::A, dna4::A, dna4::A}));

    // enlarge with value
    t0.resize(5, dna4::C);
    EXPECT_EQ(t0, (TypeParam{dna4::A, dna4::A, dna4::A, dna4::C, dna4::C}));

    // shrink with value (no effect)
    t0.resize(4, dna4::G);
    EXPECT_EQ(t0, (TypeParam{dna4::A, dna4::A, dna4::A, dna4::C}));

    // shrink without value
    t0.resize(2);
    EXPECT_EQ(t0, (TypeParam{dna4::A, dna4::A}));
}

TYPED_TEST(container, swap)
{
    TypeParam t0{};
    TypeParam t1{dna4::A, dna4::C, dna4::C, dna4::G, dna4::T};

    t0.swap(t1);
    EXPECT_EQ(t0, (TypeParam{dna4::A, dna4::C, dna4::C, dna4::G, dna4::T}));
    EXPECT_EQ(t1, (TypeParam{}));
}

TYPED_TEST(container, streamable)
{
    TypeParam t1{dna4::A, dna4::C, dna4::C, dna4::G, dna4::T};

    std::ostringstream o;
    debug_stream.set_underlying_stream(o);

    debug_stream << TypeParam{};

    o.flush();
    EXPECT_EQ(o.str(), "");

    debug_stream << t1;

    o.flush();
    EXPECT_EQ(o.str(), "ACCGT");
}

#if 0 // SEQAN3_WITH_CEREAL
template <typename in_archive_t, typename out_archive_t, typename TypeParam>
void do_serialisation(TypeParam const l)
{
    // This makes sure the file is also deleted if an exception is thrown in one of the tests below
    // Generate unique file name.
    test::tmp_file_name filename{"container_cereal_test"};
    {
        std::ofstream os{filename.get_path(), std::ios::binary};
        out_archive_t oarchive{os};
        oarchive(l);
    }

    {
        TypeParam in_l{};
        std::ifstream is{filename.get_path(), std::ios::binary};
        in_archive_t iarchive{is};
        iarchive(in_l);
        EXPECT_EQ(l, in_l);
    }
}

TYPED_TEST(container, serialisation)
{
    TypeParam t1{dna4::A, dna4::C, dna4::C, dna4::G, dna4::T};

    do_serialisation<cereal::BinaryInputArchive,         cereal::BinaryOutputArchive>        (t1);
    do_serialisation<cereal::PortableBinaryInputArchive, cereal::PortableBinaryOutputArchive>(t1);
    do_serialisation<cereal::JSONInputArchive,           cereal::JSONOutputArchive>          (t1);
    do_serialisation<cereal::XMLInputArchive,            cereal::XMLOutputArchive>           (t1);
}
#endif // SEQAN3_WITH_CEREAL

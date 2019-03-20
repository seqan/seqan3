// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>
// TODO clean up includes
#include <range/v3/algorithm/equal.hpp>
#include <range/v3/range_traits.hpp>
#include <range/v3/view/take.hpp>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/io/stream/debug_stream.hpp>
#include <seqan3/range/container/all.hpp>
#include <seqan3/range/view/convert.hpp>
#include <seqan3/search/dream_index/detail/bitvector.hpp>

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

template <typename T>
class container : public ::testing::Test
{};

using container_types = ::testing::Types<detail::bitvector<uncompressed>>;

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
    TypeParam t3{1, 0, 1, 1, 1, 0};
    TypeParam t4 = {1, 0, 1, 1, 1, 0};
    EXPECT_EQ(t3, t4);

    // n * value
    TypeParam t5{2, 1};

    // from another TypeParam's sub-range
    TypeParam t6{t3.begin() + 2, t3.begin() + 4};
    EXPECT_EQ(t5, t6);

    // direct from another container
    TypeParam t7{TypeParam{1, 0, 1, 1, 1, 0}};
    EXPECT_EQ(t3, t7);
}

TYPED_TEST(container, assign)
{
    TypeParam t0{1, 1};
    TypeParam t1{0, 1, 1, 0, 1};

    // n * value
    TypeParam t3;
    t3.assign(2, 1);
    EXPECT_EQ(t3, t0);

    // from another container's sub-range
    TypeParam t4;
    t4.assign(t1.cbegin(), t1.cend());
    EXPECT_EQ(t4, t1);

    // initializer list
    TypeParam t5, t6;
    t5.assign({0, 1, 1, 0, 1});
    t6 = {0, 1, 1, 0, 1};
    EXPECT_EQ(t5, t1);
    EXPECT_EQ(t6, t1);
}

TYPED_TEST(container, iterators)
{
    TypeParam t1{0, 1, 1, 0, 1};
    TypeParam const t2{0, 1, 1, 0, 1};

    // begin
    EXPECT_EQ(*t1.begin(),  0u);
    EXPECT_EQ(*t1.cbegin(), 0u);
    EXPECT_EQ(*t2.begin(),  0u);
    EXPECT_EQ(*t2.cbegin(), 0u);

    // end and arithmetic
    EXPECT_EQ(*(t1.end()  - 1), 1u);
    EXPECT_EQ(*(t1.cend() - 1), 1u);
    EXPECT_EQ(*(t2.end()  - 1), 1u);
    EXPECT_EQ(*(t2.cend() - 1), 1u);

    // convertibility between const and non-const
    EXPECT_TRUE(t1.cend() == t1.end());

    // mutability
    *t1.begin() = 1;
    EXPECT_TRUE((std::ranges::equal(t1, TypeParam{1, 1, 1, 0, 1})));
}

TYPED_TEST(container, element_access)
{
    TypeParam t1{0, 1, 1, 0, 1};
    TypeParam const t2{0, 1, 1, 0, 1};

    // at
    EXPECT_EQ(t1.at(0), 0u);
    EXPECT_EQ(t2.at(0), 0u);

    //TODO check at's ability to throw

    // []
    EXPECT_EQ(t1[0], 0u);
    EXPECT_EQ(t2[0], 0u);

    // front
    EXPECT_EQ(t1.front(), 0u);
    EXPECT_EQ(t2.front(), 0u);

    // back
    EXPECT_EQ(t1.back(), 1u);
    EXPECT_EQ(t2.back(), 1u);

    // mutability
    t1[0] = 1;
    EXPECT_TRUE((std::ranges::equal(t1, TypeParam{1, 1, 1, 0, 1})));

    t1.front() = 0;
    EXPECT_TRUE((std::ranges::equal(t1, TypeParam{0, 1, 1, 0, 1})));

    t1.back() = 0;
    EXPECT_TRUE((std::ranges::equal(t1, TypeParam{0, 1, 1, 0, 0})));
}

TYPED_TEST(container, capacity)
{
    TypeParam t0{};
    TypeParam t1{0, 1, 1, 0, 1};
    TypeParam const t2{0, 1, 1, 0, 1};

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
    EXPECT_LE(t1.capacity(), std::max<size_t>(t1.size()*2, 64ul)); // sdsl allocates multiples of 64
}

TYPED_TEST(container, clear)
{
    TypeParam t0{};
    TypeParam t1{0, 1, 1, 0, 1};

    t1.clear();
    EXPECT_EQ(t0, t1);
}

TYPED_TEST(container, insert)
{
    TypeParam t0{};
    TypeParam t1{0, 1, 1, 0, 1};

    // position, value
    t0.insert(t0.cend(), 0);
    t0.insert(t0.cend(), 1);
    t0.insert(t0.cend(), 0);
    t0.insert(t0.cend(), 1);
    t0.insert(t0.cbegin() + 1, 1);
    EXPECT_EQ(t0, t1);

    // position, n times values
    t0.clear();
    t0.insert(t0.cend(), 2, 1);
    t0.insert(t0.cend(), 1, 0);
    t0.insert(t0.cend(), 1, 1);
    t0.insert(t0.cbegin(), 1, 0);
    EXPECT_EQ(t0, t1);

    // iterator pair
    t0.clear();
    t0.insert(t0.cend(), t1.begin() + 1, t1.begin() + 3);

    t0.insert(t0.cend(),   t1.cend() - 2, t1.cend());
    t0.insert(t0.cbegin(), t1.cbegin(), t1.cbegin() + 1);
    EXPECT_EQ(t0, t1);

    // initializer list
    t0.clear();
    t0.insert(t0.cend(), {0, 1, 0, 1});
    t0.insert(t0.cbegin() + 1, 1);
    EXPECT_EQ(t0, t1);
}

TYPED_TEST(container, erase)
{
    TypeParam t1{0, 1, 1, 0, 1};

    // one element
    t1.erase(t1.begin());
    EXPECT_EQ(t1, (TypeParam{1, 1, 0, 1}));

    // range
    t1.erase(t1.begin() + 1, t1.begin() + 3);
    EXPECT_EQ(t1, (TypeParam{1, 1}));
}

TYPED_TEST(container, push_pop)
{
    TypeParam t0{};

    // push_back
    t0.push_back(0);
    EXPECT_EQ(t0,  (TypeParam{0}));
    t0.push_back(1);
    EXPECT_EQ(t0, (TypeParam{0, 1}));

    // pop_back
    t0.pop_back();
    EXPECT_EQ(t0, (TypeParam{0}));
    t0.pop_back();
    EXPECT_EQ(t0, (TypeParam{}));
}

TYPED_TEST(container, resize)
{
    TypeParam t0{};

    // enlarge without values
    t0.resize(3);
    EXPECT_EQ(t0, (TypeParam{0, 0, 0}));

    // enlarge with value
    t0.resize(5, 1);
    EXPECT_EQ(t0, (TypeParam{0, 0, 0, 1, 1}));

    // shrink with value (no effect)
    t0.resize(4, 0);
    EXPECT_EQ(t0, (TypeParam{0, 0, 0, 1}));

    // shrink without value
    t0.resize(2);
    EXPECT_EQ(t0, (TypeParam{0, 0}));
}

TYPED_TEST(container, swap)
{
    TypeParam t0{};
    TypeParam t1{0, 1, 1, 0, 1};

    t0.swap(t1);
    EXPECT_EQ(t0, (TypeParam{0, 1, 1, 0, 1}));
    EXPECT_EQ(t1, (TypeParam{}));
}

TYPED_TEST(container, streamable)
{
    TypeParam t1{0, 1, 1, 0, 1};

    std::ostringstream o;
    debug_stream_type my_stream{o};

    my_stream << TypeParam{};

    o.flush();
    EXPECT_EQ(o.str(), "[]");

    my_stream << t1;

    o.flush();
    EXPECT_EQ(o.str(), "[][0,1,1,0,1]");
}

/*
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
    TypeParam t1{0, 1, 1, 0, 1};

    do_serialisation<cereal::BinaryInputArchive,         cereal::BinaryOutputArchive>        (t1);
    do_serialisation<cereal::PortableBinaryInputArchive, cereal::PortableBinaryOutputArchive>(t1);
    do_serialisation<cereal::JSONInputArchive,           cereal::JSONOutputArchive>          (t1);
    do_serialisation<cereal::XMLInputArchive,            cereal::XMLOutputArchive>           (t1);
}
#endif // SEQAN3_WITH_CEREAL
*/

// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <range/v3/algorithm/equal.hpp>
#include <range/v3/range_traits.hpp>
#include <range/v3/view/take.hpp>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/container/all.hpp>
#include <seqan3/range/view/convert.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/test/cereal.hpp>
#include <seqan3/test/pretty_printing.hpp>

using namespace seqan3;

template <typename T>
class container : public ::testing::Test
{};

using container_types = ::testing::Types<std::vector<dna4>,
                                         bitcompressed_vector<dna4>,
                                         small_vector<dna4, 1000>>;

TYPED_TEST_CASE(container, container_types);

TYPED_TEST(container, concepts)
{
    EXPECT_TRUE(ReservableContainer<TypeParam>);
}

TYPED_TEST(container, construction)
{
    TypeParam t1;
    TypeParam t2{};
    EXPECT_EQ(t1, t2);

    // initializer list
    TypeParam t3{'A'_dna4, 'C'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4};
    TypeParam t4 = {'A'_dna4, 'C'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4};
    EXPECT_EQ(t3, t4);

    // n * value
    TypeParam t5{2, 'C'_dna4};

    // from another TypeParam's sub-range
    TypeParam t6{t3.begin() + 1, t3.begin() + 3};
    EXPECT_EQ(t5, t6);

    // direct from another container
    TypeParam t7{"ACCGT"_dna4};
    EXPECT_EQ(t3, t7);
}

TYPED_TEST(container, swap)
{
    TypeParam t0{};
    TypeParam t1{'A'_dna4, 'C'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4};

    t0.swap(t1);
    EXPECT_EQ(t0, (TypeParam{'A'_dna4, 'C'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4}));
    EXPECT_EQ(t1, (TypeParam{}));

    swap(t0, t1);
    EXPECT_EQ(t0, (TypeParam{}));
    EXPECT_EQ(t1, (TypeParam{'A'_dna4, 'C'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4}));
}

TYPED_TEST(container, assign)
{
    TypeParam t0{'C'_dna4, 'C'_dna4};
    TypeParam t1{'A'_dna4, 'C'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4};

    // n * value
    TypeParam t3;
    t3.assign(2, 'C'_dna4);
    EXPECT_EQ(t3, t0);

    // from another container's sub-range
    TypeParam t4;
    t4.assign(t1.cbegin(), t1.cend());
    EXPECT_EQ(t4, t1);

    // initializer list
    TypeParam t5, t6;
    t5.assign({'A'_dna4, 'C'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4});
    t6 = {'A'_dna4, 'C'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4};
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
    TypeParam t1{'A'_dna4, 'C'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4};
    TypeParam const t2{'A'_dna4, 'C'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4};

    // begin
    EXPECT_EQ(*t1.begin(),  'A'_dna4);
    EXPECT_EQ(*t1.cbegin(), 'A'_dna4);
    EXPECT_EQ(*t2.begin(),  'A'_dna4);
    EXPECT_EQ(*t2.cbegin(), 'A'_dna4);

    // end and arithmetic
    EXPECT_EQ(*(t1.end()  - 1), 'T'_dna4);
    EXPECT_EQ(*(t1.cend() - 1), 'T'_dna4);
    EXPECT_EQ(*(t2.end()  - 1), 'T'_dna4);
    EXPECT_EQ(*(t2.cend() - 1), 'T'_dna4);

    // convertibility between const and non-const
    EXPECT_TRUE(t1.cend() == t1.end());

    // mutability
    *t1.begin() = 'T'_dna4;
    EXPECT_TRUE((std::ranges::equal(t1, "TCCGT"_dna4)));
}

TYPED_TEST(container, element_access)
{
    TypeParam t1{'A'_dna4, 'C'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4};
    TypeParam const t2{'A'_dna4, 'C'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4};

    // at
    EXPECT_EQ(t1.at(0), 'A'_dna4);
    EXPECT_EQ(t2.at(0), 'A'_dna4);
    EXPECT_THROW(t1.at(20), std::out_of_range);
    EXPECT_THROW(t2.at(20), std::out_of_range);

    // []
    EXPECT_EQ(t1[0], 'A'_dna4);
    EXPECT_EQ(t2[0], 'A'_dna4);

    // front
    EXPECT_EQ(t1.front(), 'A'_dna4);
    EXPECT_EQ(t2.front(), 'A'_dna4);

    // back
    EXPECT_EQ(t1.back(), 'T'_dna4);
    EXPECT_EQ(t2.back(), 'T'_dna4);

    // mutability
    t1[0] = 'T'_dna4;
    EXPECT_TRUE((std::ranges::equal(t1, "TCCGT"_dna4)));

    t1.front() = 'C'_dna4;
    EXPECT_TRUE((std::ranges::equal(t1, "CCCGT"_dna4)));

    t1.back() = 'G'_dna4;
    EXPECT_TRUE((std::ranges::equal(t1, "CCCGG"_dna4)));
}

TYPED_TEST(container, capacity)
{
    TypeParam t0{};
    TypeParam t1{'A'_dna4, 'C'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4};
    TypeParam const t2{'A'_dna4, 'C'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4};

    // empty
    EXPECT_TRUE(t0.empty());
    EXPECT_FALSE(t1.empty());
    EXPECT_FALSE(t2.empty());

    // size
    EXPECT_EQ(t0.size(), 0u);
    EXPECT_EQ(t1.size(), 5u);
    EXPECT_EQ(t2.size(), 5u);

    // capacity
    EXPECT_GE(t0.capacity(), t0.size());
    EXPECT_GE(t1.capacity(), t1.size());
    EXPECT_GE(t2.capacity(), t2.size());

    if constexpr (!std::Same<TypeParam, small_vector<dna4, 1000>>)
    {
        // max_size
        EXPECT_GT(t0.max_size(), 1'000'000'000'000u);
        EXPECT_GT(t1.max_size(), 1'000'000'000'000u);
        EXPECT_GT(t2.max_size(), 1'000'000'000'000u);

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
    else
    {
        // max_size
        EXPECT_EQ(t0.max_size(), 1000u);
        EXPECT_EQ(t1.max_size(), 1000u);
        EXPECT_EQ(t2.max_size(), 1000u);

        // reserve
        t0.reserve(2000);
        EXPECT_EQ(t0.capacity(), 1000u); // no-op
        t1.shrink_to_fit();
        EXPECT_EQ(t0.capacity(), 1000u); // no-op
    }
}

TYPED_TEST(container, clear)
{
    TypeParam t0{};
    TypeParam t1{'A'_dna4, 'C'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4};

    t1.clear();
    EXPECT_EQ(t0, t1);
}

TYPED_TEST(container, insert)
{
    TypeParam t0{};
    TypeParam t1{'A'_dna4, 'C'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4};

    // position, value
    t0.insert(t0.cend(), 'A'_dna4);
    t0.insert(t0.cend(), 'C'_dna4);
    t0.insert(t0.cend(), 'G'_dna4);
    t0.insert(t0.cend(), 'T'_dna4);
    t0.insert(t0.cbegin() + 1, 'C'_dna4);
    EXPECT_EQ(t0, t1);

    // position, n times values
    t0.clear();
    t0.insert(t0.cend(), 2, 'C'_dna4);
    t0.insert(t0.cend(), 1, 'G'_dna4);
    t0.insert(t0.cend(), 1, 'T'_dna4);
    t0.insert(t0.cbegin(), 1, 'A'_dna4);
    EXPECT_EQ(t0, t1);

    // iterator pair
    t0.clear();
    t0.insert(t0.cend(), t1.begin() + 1, t1.begin() + 3);

    t0.insert(t0.cend(),   t1.cend() - 2, t1.cend());
    t0.insert(t0.cbegin(), t1.cbegin(), t1.cbegin() + 1);
    EXPECT_EQ(t0, t1);

    // initializer list
    t0.clear();
    t0.insert(t0.cend(), {'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4});
    t0.insert(t0.cbegin() + 1, 'C'_dna4);
    EXPECT_EQ(t0, t1);
}

TYPED_TEST(container, erase)
{
    TypeParam t1{'A'_dna4, 'C'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4};

    // one element
    t1.erase(t1.begin());
    EXPECT_EQ(t1, (TypeParam{'C'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4}));

    // range
    t1.erase(t1.begin() + 1, t1.begin() + 3);
    EXPECT_EQ(t1, (TypeParam{'C'_dna4, 'T'_dna4}));

    // empty range (no op)
    t1.erase(t1.begin(), t1.begin());
    EXPECT_EQ(t1, (TypeParam{'C'_dna4, 'T'_dna4}));
}

TYPED_TEST(container, push_pop)
{
    TypeParam t0{};

    // push_back
    t0.push_back('A'_dna4);
    EXPECT_EQ(t0,  (TypeParam{'A'_dna4}));
    t0.push_back('C'_dna4);
    EXPECT_EQ(t0, (TypeParam{'A'_dna4, 'C'_dna4}));

    // pop_back
    t0.pop_back();
    EXPECT_EQ(t0, (TypeParam{'A'_dna4}));
    t0.pop_back();
    EXPECT_EQ(t0, (TypeParam{}));
}

TYPED_TEST(container, resize)
{
    TypeParam t0{};

    // enlarge without values
    t0.resize(3);
    EXPECT_EQ(t0, (TypeParam{'A'_dna4, 'A'_dna4, 'A'_dna4}));

    // enlarge with value
    t0.resize(5, 'C'_dna4);
    EXPECT_EQ(t0, (TypeParam{'A'_dna4, 'A'_dna4, 'A'_dna4, 'C'_dna4, 'C'_dna4}));

    // shrink with value (no effect)
    t0.resize(4, 'G'_dna4);
    EXPECT_EQ(t0, (TypeParam{'A'_dna4, 'A'_dna4, 'A'_dna4, 'C'_dna4}));

    // shrink without value
    t0.resize(2);
    EXPECT_EQ(t0, (TypeParam{'A'_dna4, 'A'_dna4}));
}

TYPED_TEST(container, streamable)
{
    TypeParam t1{'A'_dna4, 'C'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4};

    std::ostringstream o;
    debug_stream_type my_stream{o};

    my_stream << TypeParam{};

    o.flush();
    EXPECT_EQ(o.str(), "");

    my_stream << t1;

    o.flush();
    EXPECT_EQ(o.str(), "ACCGT");
}

TYPED_TEST(container, serialisation)
{
    TypeParam t1{'A'_dna4, 'C'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4};
    test::do_serialisation(t1);
}

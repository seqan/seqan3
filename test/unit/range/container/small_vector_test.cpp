// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <string>
#include <utility>

#include <seqan3/range/container/concept.hpp>
#include <seqan3/range/container/small_vector.hpp>
#include <seqan3/std/ranges>
#include <seqan3/test/cereal.hpp>

// standard construction.
TEST(small_vector, standard_construction)
{
    EXPECT_TRUE((std::is_default_constructible_v<seqan3::small_vector<char, 4>>));
    EXPECT_TRUE((std::is_nothrow_default_constructible_v<seqan3::small_vector<char, 4>>));
    EXPECT_TRUE((std::is_copy_constructible_v<seqan3::small_vector<char, 4>>));
    EXPECT_TRUE((std::is_trivially_copy_constructible_v<seqan3::small_vector<char, 4>>));
    EXPECT_TRUE((std::is_nothrow_copy_constructible_v<seqan3::small_vector<char, 4>>));
    EXPECT_TRUE((std::is_move_constructible_v<seqan3::small_vector<char, 4>>));
    EXPECT_TRUE((std::is_trivially_move_constructible_v<seqan3::small_vector<char, 4>>));
    EXPECT_TRUE((std::is_nothrow_move_constructible_v<seqan3::small_vector<char, 4>>));
    EXPECT_TRUE((std::is_copy_assignable_v<seqan3::small_vector<char, 4>>));
    EXPECT_TRUE((std::is_trivially_copy_assignable_v<seqan3::small_vector<char, 4>>));
    EXPECT_TRUE((std::is_nothrow_copy_assignable_v<seqan3::small_vector<char, 4>>));
    EXPECT_TRUE((std::is_move_assignable_v<seqan3::small_vector<char, 4>>));
    EXPECT_TRUE((std::is_trivially_move_assignable_v<seqan3::small_vector<char, 4>>));
    EXPECT_TRUE((std::is_nothrow_move_assignable_v<seqan3::small_vector<char, 4>>));
}

TEST(small_vector, concepts)
{
    EXPECT_TRUE((seqan3::reservible_container<seqan3::small_vector<char, 4>>));
    EXPECT_TRUE((std::ranges::random_access_range<seqan3::small_vector<char, 4>>));
    EXPECT_TRUE((std::ranges::contiguous_range<seqan3::small_vector<char, 4>>));
}

TEST(small_vector, construct_from_array)
{
    // Deduce value type and N
    EXPECT_TRUE((std::same_as<decltype(seqan3::small_vector{std::array{'h','e','l','l','o'}}),
                 seqan3::small_vector<char, 5>>));

    // construct from different sized array (size has to be smaller)
    EXPECT_TRUE((std::same_as<decltype(seqan3::small_vector<char, 10>{std::array{'h','e','l','l','o'}}),
                 seqan3::small_vector<char, 10>>));
}

TEST(small_vector, construct_from_built_in_array)
{
    // Deduce value type and N
    int arr[3] = {1, 2, 3};
    EXPECT_TRUE((std::same_as<decltype(seqan3::small_vector{arr}), seqan3::small_vector<int, 3>>));

    // Deduce value type and N
    EXPECT_TRUE((std::same_as<decltype(seqan3::small_vector<int, 5>{arr}), seqan3::small_vector<int, 5>>));

    // char const *
    EXPECT_TRUE((std::same_as<decltype(seqan3::small_vector{"hi"}), seqan3::small_vector<char, 3>>));

    // parameter pack
    EXPECT_TRUE((std::same_as<decltype(seqan3::small_vector<char, 3>{'A', 'C', 'X'}), seqan3::small_vector<char, 3>>));
}

// ---------------------------------------------------------------------------------------------------------------------
// constexpr tests
// ---------------------------------------------------------------------------------------------------------------------
constexpr bool comparison_test()
{
    seqan3::small_vector<char, 20> t1{'A', 'C', 'C', 'G', 'T'};
    seqan3::small_vector<char, 20> t2{'A', 'C', 'C', 'G', 'T'};
    seqan3::small_vector<char, 2> t3{'A', 'C'};
    seqan3::small_vector<char, 20> t4{'A', 'G', 'C', 'G', 'T'};

    bool res = t1 == t2;
    res = res && (t1 <= t2);
    res = res && (t1 >= t2);
    res = res && (t1 != t3);

    res = res && (t3 <  t1);
    res = res && (t3 <= t1);
    res = res && (t1 <  t4);
    res = res && (t1 <= t4);

    res = res && (t1 >  t3);
    res = res && (t1 >= t3);
    res = res && (t4 >  t1);
    res = res && (t4 >= t1);

    return true;
}


TEST(small_vector, comparison)
{
    constexpr bool b = comparison_test();
    EXPECT_TRUE(b);
}

constexpr bool begin_end_test()
{
    std::array src{'h','e','l','l','o'};
    seqan3::small_vector vec{src};

    auto it_s = src.begin();
    auto it_v = vec.begin();

    for (; it_v != vec.end(); ++it_s, ++it_v)
        if (*it_v != *it_s)
            return false;

    return true;
}

constexpr bool cbegin_cend_test()
{
    std::array src{'h','e','l','l','o'};
    seqan3::small_vector vec{src};

    auto it_s = src.cbegin();
    auto it_v = vec.cbegin();

    for (; it_v != vec.cend(); ++it_s, ++it_v)
        if (*it_v != *it_s)
            return false;

    return true;
}

TEST(small_vector, iterator)
{
    constexpr bool b_const = begin_end_test();
    EXPECT_TRUE(b_const);

    constexpr bool cb_const = cbegin_cend_test();
    EXPECT_TRUE(cb_const);
}

TEST(small_string, size_and_maxsize)
{
    // auto deduction -> capacity == size
    {
        constexpr seqan3::small_vector vec{"hello"};
        constexpr auto size = vec.size();
        constexpr auto msize = vec.max_size();

        EXPECT_EQ(size, 6u); // incl. null character
        EXPECT_EQ(msize, 6u);
    }

    // capacity != size
    {
        constexpr seqan3::small_vector<char, 10> vec{'h','e','l','l','o'};
        constexpr auto size = vec.size();
        constexpr auto msize = vec.max_size();

        EXPECT_EQ(size, 5u);
        EXPECT_EQ(msize, 10u);
    }
}

constexpr bool swap_test()
{
    seqan3::small_vector<char, 20> t0{};
    seqan3::small_vector<char, 20> t1{"AC"};

    t0.swap(t1);
    return (t0 == seqan3::small_vector<char, 20>{"AC"}) && (t1 == seqan3::small_vector<char, 20>{});
}

TEST(small_vector, swap)
{
    constexpr bool res = swap_test();
    EXPECT_TRUE(res);
}

constexpr bool assign_test()
{
    seqan3::small_vector<char, 20> t0{'C', 'C'};
    seqan3::small_vector<char, 20> t1{'A', 'C', 'C', 'G', 'T'};

    // n * value
    seqan3::small_vector<char, 20> t3;
    t3.assign(2, 'C');
    bool res = (t3 == t0);

    // from another container's sub-range
    seqan3::small_vector<char, 20> t4;
    t4.assign(t1.cbegin(), t1.cend());
    res = res && (t4 == t1);

    // initializer list
    seqan3::small_vector<char, 20> t5, t6;
    t5.assign({'A', 'C', 'C', 'G', 'T'});
    t6 = {'A', 'C', 'C', 'G', 'T'};
    res = res && (t5 == t1);
    res = res && (t6 == t1);

    // direct from another container
    seqan3::small_vector<char, 20> t7;
    t7.assign(std::array<char, 5>{'A', 'C', 'C', 'G', 'T'});
    res = res && (t7 == t1);

    return res;
}

TEST(small_vector, assign)
{
    constexpr bool res = assign_test();
    EXPECT_TRUE(res);
}

constexpr bool element_access_test()
{
    seqan3::small_vector<char, 20> t1{'A', 'C', 'C', 'G', 'T'};
    seqan3::small_vector<char, 20> const t2{'A', 'C', 'C', 'G', 'T'};

    // at() cannot be constexpr because it throws

    // []
    bool res = (t1[0] == 'A');
    res = res && (t2[0] == 'A');
    // front
    res = res && (t1.front() == 'A');
    res = res && (t2.front() == 'A');

    // back
    res = res && (t1.back() == 'T');
    res = res && (t2.back() == 'T');

    // mutability
    t1[0] = 'T';
    res = res && (t1 == seqan3::small_vector<char, 20>{'T', 'C', 'C', 'G', 'T'});
    res = res && !(t1 == seqan3::small_vector<char, 20>{'T', 'C', 'C'}); // for code coverage

    t1.front() = 'C';
    res = res && (t1 == seqan3::small_vector<char, 20>{'C', 'C', 'C', 'G', 'T'});

    t1.back() = 'G';
    res = res && (t1 == seqan3::small_vector<char, 20>{'C', 'C', 'C', 'G', 'G'});

    // data()
    res = res && (*t1.data() == 'C');
    res = res && (*t2.data() == 'A');

    return res;
}

TEST(small_vector, element_access)
{
    constexpr bool res = element_access_test();
    EXPECT_TRUE(res);
}

constexpr bool clear_test()
{
    seqan3::small_vector<char, 20> t0{};
    seqan3::small_vector<char, 20> t1{'A', 'C', 'C', 'G', 'T'};

    t1.clear();

    return t0 == t1;
}

TEST(small_vector, clear)
{
    constexpr bool res = clear_test();
    EXPECT_TRUE(res);
}

constexpr bool insert_test()
{
    seqan3::small_vector<char, 20> t0{};
    seqan3::small_vector<char, 20> t1{'A', 'C', 'C', 'G', 'T'};

    // position, value
    t0.insert(t0.cend(), 'A');
    t0.insert(t0.cend(), 'C');
    t0.insert(t0.cend(), 'G');
    t0.insert(t0.cend(), 'T');
    t0.insert(t0.cbegin() + 1, 'C');
    bool res = (t0 == t1);

    // position, n times values
    t0.clear();
    t0.insert(t0.cend(), 2, 'C');
    t0.insert(t0.cend(), 1, 'G');
    t0.insert(t0.cend(), 1, 'T');
    t0.insert(t0.cbegin(), 1, 'A');
    res = res && (t0 == t1);

    // iterator pair
    t0.clear();
    t0.insert(t0.cend(), t1.begin() + 1, t1.begin() + 3);

    t0.insert(t0.cend(),   t1.cend() - 2, t1.cend());
    t0.insert(t0.cbegin(), t1.cbegin(), t1.cbegin() + 1);
    res = res && (t0 == t1);

    // initializer list
    t0.clear();
    t0.insert(t0.cend(), {'A', 'C', 'G', 'T'});
    t0.insert(t0.cbegin() + 1, 'C');
    res = res && (t0 == t1);

    return res;
}

TEST(small_vector, insert)
{
    constexpr bool res = insert_test();
    EXPECT_TRUE(res);
}

constexpr bool erase_test()
{
    seqan3::small_vector<char, 20> t1{'A', 'C', 'C', 'G', 'T'};

    // one element
    t1.erase(t1.begin());
    bool res = (t1 == (seqan3::small_vector<char, 20>{'C', 'C', 'G', 'T'}));

    // range
    t1.erase(t1.begin() + 1, t1.begin() + 3);
    res = res && (t1 == (seqan3::small_vector<char, 20>{'C', 'T'}));

    return res;
}

TEST(small_vector, erase)
{
    constexpr bool res = erase_test();
    EXPECT_TRUE(res);
}

constexpr bool push_pop_test()
{
    seqan3::small_vector<char, 20> t0{};

    // push_back
    t0.push_back('A');
    bool res = (t0 ==  (seqan3::small_vector<char, 20>{'A'}));
    t0.push_back('C');
    res = res && (t0 == (seqan3::small_vector<char, 20>{'A', 'C'}));

    // pop_back
    t0.pop_back();
    res = res && (t0 == (seqan3::small_vector<char, 20>{'A'}));
    t0.pop_back();
    res = res && (t0 == (seqan3::small_vector<char, 20>{}));

    return res;
}

TEST(small_vector, push_pop)
{
    constexpr bool res = push_pop_test();
    EXPECT_TRUE(res);
}

constexpr bool resize_test()
{
    seqan3::small_vector<int, 20> t0{};

    // enlarge without values
    t0.resize(3);
    bool res = (t0 == (seqan3::small_vector<int, 20>{0, 0, 0}));

    // enlarge with value
    t0.resize(5, 11);
    res = res && (t0 == (seqan3::small_vector<int, 20>{0, 0, 0, 11, 11}));

    // shrink with value (no effect)
    t0.resize(4, 500);
    res = res && (t0 == (seqan3::small_vector<int, 20>{0, 0, 0, 11}));

    // shrink without value
    t0.resize(2);
    res = res && (t0 == (seqan3::small_vector<int, 20>{0, 0}));

    return res;
}

TEST(small_vector, resize)
{
    constexpr bool res = resize_test();
    EXPECT_TRUE(res);
}

TEST(small_vector, serialisation)
{
    seqan3::small_vector hello{std::array{'h','e','l','l','o'}};
    seqan3::test::do_serialisation(hello);
}

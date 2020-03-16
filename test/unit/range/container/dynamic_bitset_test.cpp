// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/container/dynamic_bitset.hpp>
#include <seqan3/test/cereal.hpp>

// Standard construction.
TEST(dynamic_bitset, standard_construction)
{
    EXPECT_TRUE((std::is_default_constructible_v<seqan3::dynamic_bitset<58>>));
    EXPECT_TRUE((std::is_nothrow_default_constructible_v<seqan3::dynamic_bitset<58>>));
    EXPECT_TRUE((std::is_copy_constructible_v<seqan3::dynamic_bitset<58>>));
    EXPECT_TRUE((std::is_trivially_copy_constructible_v<seqan3::dynamic_bitset<58>>));
    EXPECT_TRUE((std::is_nothrow_copy_constructible_v<seqan3::dynamic_bitset<58>>));
    EXPECT_TRUE((std::is_move_constructible_v<seqan3::dynamic_bitset<58>>));
    EXPECT_TRUE((std::is_trivially_move_constructible_v<seqan3::dynamic_bitset<58>>));
    EXPECT_TRUE((std::is_nothrow_move_constructible_v<seqan3::dynamic_bitset<58>>));
    EXPECT_TRUE((std::is_copy_assignable_v<seqan3::dynamic_bitset<58>>));
    EXPECT_TRUE((std::is_trivially_copy_assignable_v<seqan3::dynamic_bitset<58>>));
    EXPECT_TRUE((std::is_nothrow_copy_assignable_v<seqan3::dynamic_bitset<58>>));
    EXPECT_TRUE((std::is_move_assignable_v<seqan3::dynamic_bitset<58>>));
    EXPECT_TRUE((std::is_trivially_move_assignable_v<seqan3::dynamic_bitset<58>>));
    EXPECT_TRUE((std::is_nothrow_move_assignable_v<seqan3::dynamic_bitset<58>>));
    EXPECT_THROW(seqan3::dynamic_bitset{std::numeric_limits<uint64_t>::max()}, std::invalid_argument);
    EXPECT_THROW(seqan3::dynamic_bitset{"10101011x0101"}, std::invalid_argument);
    EXPECT_EQ(seqan3::detail::sizeof_bits<decltype(*std::declval<seqan3::dynamic_bitset<58>>().raw_data())>, 64u);
}

TEST(dynamic_bitset, concepts)
{
    EXPECT_TRUE((seqan3::reservible_container<seqan3::dynamic_bitset<58>>));
    EXPECT_TRUE((std::ranges::random_access_range<seqan3::dynamic_bitset<58>>));
}

constexpr bool comparison_test()
{
    constexpr seqan3::dynamic_bitset t1{0b11'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111};
    constexpr seqan3::dynamic_bitset t2{0b11'1111'1111'1111'1111'1111'1100'1111'1111'1111'1111'1111'1111'1111'1100};
    constexpr seqan3::dynamic_bitset t3{0b00'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111};
    constexpr seqan3::dynamic_bitset t4{72'057'594'037'927'935};
    constexpr seqan3::dynamic_bitset t5{"1111111111111111111111111111111111111111111111111111111111"};
    constexpr seqan3::dynamic_bitset t6{"1111111111111111111111110011111111111111111111111111111100"};
    constexpr seqan3::dynamic_bitset t7{"11111111111111111111111111111111111111111111111111111111"};

    bool res = t3 == t4;
    res &= t1 == t5;
    res &= t2 == t6;
    res &= t3 == t7;

    res &= t1 > t2;
    res &= t2 > t3;
    res &= t1 > t3;

    res &= t1 >= t2;
    res &= t2 >= t3;
    res &= t1 >= t3;

    res &= t2 <= t1;
    res &= t3 <= t2;
    res &= t3 <= t1;

    res &= t2 < t1;
    res &= t3 < t2;
    res &= t3 < t1;

    res &= t1 != t2;
    res &= t2 != t3;
    res &= t1 != t3;

    return res;
}

TEST(dynamic_bitset, comparison)
{
    constexpr bool b = comparison_test();
    EXPECT_TRUE(b);
    EXPECT_TRUE(comparison_test());
}

constexpr bool size_test()
{
    constexpr seqan3::dynamic_bitset t1{0b11'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111};
    bool res = t1.size() == 58u;

    constexpr seqan3::dynamic_bitset t2{0b11'1111'1111'1111'1111'1111'1100'1111'1111'1111'1111'1111'1111'1111'1111};
    res &= t2.size() == 58u;

    constexpr seqan3::dynamic_bitset t3{0b00'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111};
    res &= t3.size() == 56u;

    constexpr seqan3::dynamic_bitset t4{0b11'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1100};
    res &= t4.size() == 58u;

    constexpr seqan3::dynamic_bitset t5;
    res &= t5.size() == 0u;

    return res;
}

TEST(dynamic_bitset, size)
{
    constexpr bool b = size_test();
    EXPECT_TRUE(b);
    EXPECT_TRUE(size_test());
}

constexpr bool count_test()
{
    constexpr seqan3::dynamic_bitset t1{0b11'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111};
    bool res = t1.count() == 58u;

    constexpr seqan3::dynamic_bitset t2{0b11'1111'1111'1111'1111'1111'1100'1111'1111'1111'1111'1111'1111'1111'1111};
    res &= t2.count() == 56u;

    constexpr seqan3::dynamic_bitset t3{0b00'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111};
    res &= t3.count() == 56u;

    constexpr seqan3::dynamic_bitset t4{0b11'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1100};
    res &= t4.count() == 56u;

    constexpr seqan3::dynamic_bitset t5;
    res &= t5.count() == 0u;

    return res;
}

TEST(dynamic_bitset, count)
{
    constexpr bool b = count_test();
    EXPECT_TRUE(b);
    EXPECT_TRUE(count_test());
}

constexpr bool all_test()
{
    constexpr seqan3::dynamic_bitset t1{0b11'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111};
    bool res = t1.all();

    constexpr seqan3::dynamic_bitset t2{0b11'1111'1111'1111'1111'1111'1100'1111'1111'1111'1111'1111'1111'1111'1111};
    res &= !t2.all();

    constexpr seqan3::dynamic_bitset t3{0b00'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111};
    res &= t3.all();

    constexpr seqan3::dynamic_bitset t4{0b11'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1100};
    res &= !t4.all();

    constexpr seqan3::dynamic_bitset t5;
    res &= t5.all();

    return res;
}

TEST(dynamic_bitset, all)
{
    constexpr bool b = all_test();
    EXPECT_TRUE(b);
    EXPECT_TRUE(all_test());
}

constexpr bool any_test()
{
    constexpr seqan3::dynamic_bitset t1{0b11'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111};
    bool res = t1.any();

    constexpr seqan3::dynamic_bitset t2{0b11'1111'1111'1111'1111'1111'1100'1111'1111'1111'1111'1111'1111'1111'1111};
    res &= t2.any();

    constexpr seqan3::dynamic_bitset t3{0b00'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111};
    res &= t3.any();

    constexpr seqan3::dynamic_bitset t4{0b11'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1100};
    res &= t4.any();

    constexpr seqan3::dynamic_bitset t5;
    res &= !t5.any();

    return res;
}

TEST(dynamic_bitset, any)
{
    constexpr bool b = any_test();
    EXPECT_TRUE(b);
    EXPECT_TRUE(any_test());
}

constexpr bool none_test()
{
    constexpr seqan3::dynamic_bitset t1{0b11'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111};
    bool res = !t1.none();

    constexpr seqan3::dynamic_bitset t2{0b11'1111'1111'1111'1111'1111'1100'1111'1111'1111'1111'1111'1111'1111'1111};
    res &= !t2.none();

    constexpr seqan3::dynamic_bitset t3{0b00'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111};
    res &= !t3.none();

    constexpr seqan3::dynamic_bitset t4{0b11'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1100};
    res &= !t4.none();

    constexpr seqan3::dynamic_bitset t5;
    res &= t5.none();

    return res;
}

TEST(dynamic_bitset, none)
{
    constexpr bool b = none_test();
    EXPECT_TRUE(b);
    EXPECT_TRUE(none_test());
}

constexpr bool set_test()
{
    seqan3::dynamic_bitset t1{0b11'1110'1111'1111'1111'1111'1100'1111'1111'1111'1111'1111'1111'1111'1110};
    bool res = !t1.all();
    t1.set();
    res &= t1.all();

    seqan3::dynamic_bitset t2{0b1110'1111'1111'1111'1111'1100'1111'1111'1111'1111'1111'1111'1111'1110};
    res &= !t2.all();
    t2.set();
    res &= t2.all();

    seqan3::dynamic_bitset t3{0b11'1111'1111'1111'1111'1111'1101'1111'1111'1111'1111'1111'1111'1111'1110};
    res &= !t3.all();
    t3.set(0, true);
    t3.set(33, true);
    res &= t3.all();

    return res;
}

TEST(dynamic_bitset, set)
{
    constexpr bool b = set_test();
    EXPECT_TRUE(b);
    EXPECT_TRUE(set_test());
}

constexpr bool reset_test()
{
    seqan3::dynamic_bitset t1{0b11'1110'1111'1111'1111'1111'1100'1111'1111'1111'1111'1111'1111'1111'1110};
    bool res = !t1.none();
    t1.reset();
    res &= t1.none();

    seqan3::dynamic_bitset t2{0b10'0000'0000'0000'0000'0000'0000'0000'0000'0000'0000'0000'0000'0000'0001};
    res &= !t2.none();
    t2.reset(0);
    t2.reset(57);
    res &= t2.none();

    return res;
}

TEST(dynamic_bitset, reset)
{
    constexpr bool b = reset_test();
    EXPECT_TRUE(b);
    EXPECT_TRUE(reset_test());
}

constexpr bool flip_test()
{
    seqan3::dynamic_bitset t1{0b1111111111111111111111111111111111111111111111111111111111};
    bool res = t1.all();
    t1.flip();
    res &= t1.none();

    seqan3::dynamic_bitset t2{0b1111111111111111111111111111111111111111111111111111111111};
    res &= t2.all();
    t2.flip(0);
    res &= !t2.all();
    res &= t2.any();

    return res;
}

TEST(dynamic_bitset, flip)
{
    constexpr bool b = flip_test();
    EXPECT_TRUE(b);
    EXPECT_TRUE(flip_test());
}

constexpr bool access_test()
{
    bool res = true;

    seqan3::dynamic_bitset t1{0b1111'0000'0000'0000};
    seqan3::dynamic_bitset const t2{0b1111'0000'0000'0000};
    seqan3::dynamic_bitset expected{"0111000000000001"};

    for (size_t i = 0u; i < t1.size() - 4u; ++i)
    {
        res &= t1.at(i) == 0 && t1.test(i) == 0 && t1[i] == 0;
        res &= t2.at(i) == 0 && t2.test(i) == 0 && t2[i] == 0;
    }

    for (size_t i = t1.size() - 4u; i < t1.size(); ++i)
    {
        res &= t1.at(i) == 1 && t1.test(i) == 1 && t1[i] == 1;
        res &= t2.at(i) == 1 && t2.test(i) == 1 && t2[i] == 1;
    }

    res &= !t1.front();
    res &= t1.back();
    res &= !t2.front();
    res &= t2.back();

    t1[1] = true;
    res &= t1 == seqan3::dynamic_bitset{0b1111'0000'0000'0010};
    t1.at(1) = false;
    res &= t1 == t2;

    t1.front() = true;
    res &= t1 == seqan3::dynamic_bitset{0b1111'0000'0000'0001};

    t1.back() = false;
    res &= t1 == expected;

    return res;
}

TEST(dynamic_bitset, access)
{
    constexpr bool b = access_test();
    EXPECT_TRUE(b);
    EXPECT_TRUE(access_test());

    constexpr seqan3::dynamic_bitset t1{0b1111'0000'0000'0000};
    EXPECT_THROW(t1.at(16), std::out_of_range);
    EXPECT_THROW(t1.test(16), std::out_of_range);
}

constexpr bool bitwise_and_test()
{
    seqan3::dynamic_bitset t1{0b1111'0001'0000'1100};
    seqan3::dynamic_bitset t2{0b1010'0001'0000'0011};
    seqan3::dynamic_bitset expected{0b1010'0001'0000'0000};

    bool res = (t1 & t2) == expected;
    t1 &= t2;
    res = t1 == expected;

    return res;
}

TEST(dynamic_bitset, bitwise_and)
{
    constexpr bool b = bitwise_and_test();
    EXPECT_TRUE(b);
    EXPECT_TRUE(bitwise_and_test());
}

constexpr bool bitwise_or_test()
{
    seqan3::dynamic_bitset t1{0b1111'0001'0000'1100};
    seqan3::dynamic_bitset t2{0b1010'0001'0000'0011};
    seqan3::dynamic_bitset expected{0b1111'0001'0000'1111};

    bool res = (t1 | t2) == expected;
    t1 |= t2;
    res = t1 == expected;

    return res;
}

TEST(dynamic_bitset, bitwise_or)
{
    constexpr bool b = bitwise_or_test();
    EXPECT_TRUE(b);
    EXPECT_TRUE(bitwise_or_test());
}

constexpr bool bitwise_xor_test()
{
    seqan3::dynamic_bitset t1{0b1111'0001'0000'1100};
    seqan3::dynamic_bitset t2{0b1010'0001'0000'0011};
    seqan3::dynamic_bitset expected{"0101000000001111"};

    bool res = (t1 ^ t2) == expected;
    t1 ^= t2;
    res = t1 == expected;

    return res;
}

TEST(dynamic_bitset, bitwise_xor)
{
    constexpr bool b = bitwise_xor_test();
    EXPECT_TRUE(b);
    EXPECT_TRUE(bitwise_xor_test());
}

constexpr bool bitwise_not_test()
{
    seqan3::dynamic_bitset t1{0b1111'0001'0000'1100};
    seqan3::dynamic_bitset expected{"0000111011110011"};

    return ~t1 == expected;
}

TEST(dynamic_bitset, bitwise_not)
{
    constexpr bool b = bitwise_not_test();
    EXPECT_TRUE(b);
    EXPECT_TRUE(bitwise_not_test());
}

constexpr bool shift_left_test()
{
    seqan3::dynamic_bitset t1{0b1111'0001'0000'1100};

    bool res = t1 << 3 == seqan3::dynamic_bitset{0b1000'1000'0110'0000};
    t1 <<= 4;
    res &= t1 == seqan3::dynamic_bitset{"0001000011000000"};

    return res;
}

TEST(dynamic_bitset, shift_left)
{
    constexpr bool b = shift_left_test();
    EXPECT_TRUE(b);
    EXPECT_TRUE(shift_left_test());
}

constexpr bool shift_right_test()
{
    seqan3::dynamic_bitset t1{0b1111'0001'0000'1100};

    bool res = t1 >> 3 == seqan3::dynamic_bitset{"0001111000100001"};
    t1 >>= 4;
    res &= t1 == seqan3::dynamic_bitset{"0000111100010000"};

    return res;
}

TEST(dynamic_bitset, shift_right)
{
    constexpr bool b = shift_right_test();
    EXPECT_TRUE(b);
    EXPECT_TRUE(shift_right_test());
}

constexpr bool swap_test()
{
    seqan3::dynamic_bitset t1{0b1111'0001'0000'1100};
    seqan3::dynamic_bitset t2;
    seqan3::dynamic_bitset expected{t1};

    t1.swap(t2);
    bool res = t2 == expected;

    swap(t1, t2);
    res &= t1 == expected;

    t2.swap(std::move(t1));
    res &= t2 == expected;

    return res;
}

TEST(dynamic_bitset, swap)
{
    constexpr bool b = swap_test();
    EXPECT_TRUE(b);
    EXPECT_TRUE(swap_test());
}

constexpr bool assign_test()
{
    seqan3::dynamic_bitset t1{0b1111};
    seqan3::dynamic_bitset t2{0b1001};
    seqan3::dynamic_bitset t3, t4, t5, t6;

    t3.assign(4, true);
    bool res = t3 == t1;

    t4.assign(t2.cbegin(), t2.cend());
    res &= t4 == t2;

    t5.assign({true, true, true, true});
    res &= t5 == t1;

    t6 = {1, 0, 0, 1};
    res &= t6 == t2;

    return res;
}

TEST(dynamic_bitset, assign)
{
    constexpr bool b = assign_test();
    EXPECT_TRUE(b);
    EXPECT_TRUE(assign_test());
}

constexpr bool iterator_test()
{
    seqan3::dynamic_bitset t1{0b1111'0001'0000'1100};
    seqan3::dynamic_bitset const t2{0b1010'0001'0000'0011};

    bool res = !*t1.begin();
    res &= !*t1.cbegin();
    res &= *t2.begin();
    res &= *t2.cbegin();

    res &= *(t1.end() - 1);
    res &= *(t1.cend() - 1);
    res &= *(t2.end() - 1);
    res &= *(t2.cend() - 1);

    res &= t1.end() == t1.cend();

    *t1.begin() = true;
    res &= t1 == seqan3::dynamic_bitset{0b1111'0001'0000'1101};

    return res;
}

TEST(dynamic_bitset, iterators)
{
    constexpr bool b = iterator_test();
    EXPECT_TRUE(b);
    EXPECT_TRUE(iterator_test());
}

constexpr bool capacity_test()
{
    seqan3::dynamic_bitset t0{};
    seqan3::dynamic_bitset t1{0b1111'0001'0000'1100};
    seqan3::dynamic_bitset const t2{0b1010'0001'0000'0011};

    bool res = t0.empty();
    res &= !t1.empty();
    res &= !t2.empty();

    res &= t0.size() == 0u;
    res &= t1.size() == 16u;
    res &= t2.size() == 16u;

    res &= t0.capacity() >= t0.size();
    res &= t1.capacity() >= t1.size();
    res &= t2.capacity() >= t2.size();

    res &= t0.max_size() == t0.capacity();
    res &= t1.max_size() == t1.capacity();
    res &= t2.max_size() == t2.capacity();

    size_t cap = t0.capacity();
    t0.reserve(1000);
    res &= t0.capacity() == cap;
    t0.shrink_to_fit();
    res &= t0.capacity() == cap;

    res &= t1.capacity() == 58u;
    seqan3::dynamic_bitset<30> t3;
    res &= t3.capacity() == 30u;

    return res;
}

TEST(dynamic_bitset, capacity)
{
    constexpr bool b = capacity_test();
    EXPECT_TRUE(b);
    EXPECT_TRUE(capacity_test());
}

constexpr bool clear_test()
{
    seqan3::dynamic_bitset t1{0b1111'0001'0000'1100};

    t1.clear();

    return t1 == seqan3::dynamic_bitset{};
}

TEST(dynamic_bitset, clear)
{
    constexpr bool b = clear_test();
    EXPECT_TRUE(b);
    EXPECT_TRUE(clear_test());
}

constexpr bool insert_test()
{
    seqan3::dynamic_bitset t0{};
    seqan3::dynamic_bitset t1{0b100101};

    t0.insert(t0.cend(), 1);
    t0.insert(t0.cend(), 0);
    t0.insert(t0.cend(), 1);
    t0.insert(t0.cend(), 0);
    t0.insert(t0.cend(), 1);
    t0.insert(t0.cbegin() + 3, 0);
    bool res = t0 == t1;

    t0.clear();
    t0.insert(t0.cend(), 3, true);
    t0.insert(t0.cbegin() + 1, false);
    t0.insert(t0.cbegin() + 3, 2, false);
    t0.insert(t0.cbegin() + 3, 0, false);
    res &= t0 == t1;

    return res;
}

TEST(dynamic_bitset, insert)
{
    constexpr bool b = insert_test();
    EXPECT_TRUE(b);
    EXPECT_TRUE(insert_test());
}

constexpr bool erase_test()
{
    seqan3::dynamic_bitset t1{0b100101};

    t1.erase(t1.begin());
    bool res = t1 == seqan3::dynamic_bitset{0b10010};

    t1.erase(t1.begin() + 1, t1.begin() + 3);
    res &= t1 == seqan3::dynamic_bitset{0b100};

    t1.erase(t1.begin(), t1.begin());
    res &= t1 == seqan3::dynamic_bitset{0b100};

    return res;
}

TEST(dynamic_bitset, erase)
{
    constexpr bool b = erase_test();
    EXPECT_TRUE(b);
    EXPECT_TRUE(erase_test());
}

constexpr bool push_pop_test()
{
    seqan3::dynamic_bitset t1{};
    seqan3::dynamic_bitset expected{0b01};
    expected.resize(2);

    t1.push_back(true);
    bool res = t1 == seqan3::dynamic_bitset{0b1};

    t1.push_back(false);
    res &= t1 == expected;

    t1.pop_back();
    res &= t1 == seqan3::dynamic_bitset{0b1};

    t1.pop_back();
    res &= t1 == seqan3::dynamic_bitset{};

    return res;
}

TEST(dynamic_bitset, push_pop)
{
    constexpr bool b = push_pop_test();
    EXPECT_TRUE(b);
    EXPECT_TRUE(push_pop_test());
}

constexpr bool resize_test()
{
    seqan3::dynamic_bitset t1{};

    t1.resize(2);
    bool res = !t1.at(0) && !t1.at(1);

    t1.resize(5, true);
    res &= t1 == seqan3::dynamic_bitset{0b11100};

    t1.resize(4, true);
    res &= t1 == seqan3::dynamic_bitset{0b1100};

    t1.resize(3);
    res &= t1 == seqan3::dynamic_bitset{0b100};

    return res;
}

TEST(dynamic_bitset, resize)
{
    constexpr bool b = resize_test();
    EXPECT_TRUE(b);
    EXPECT_TRUE(resize_test());
}

TEST(dynamic_bitset, to_string)
{
    seqan3::dynamic_bitset t1{"0011000"};
    EXPECT_EQ(t1.to_string(), std::string{"0011000"});

    seqan3::dynamic_bitset t2{0b001100};
    EXPECT_EQ(t2.to_string(), std::string{"1100"});
    t2.resize(6);
    EXPECT_EQ(t2.to_string(), std::string{"001100"});
    EXPECT_EQ(t2.to_string('#'), std::string{"##11##"});
    EXPECT_EQ(t2.to_string('#', '*'), std::string{"##**##"});
}

constexpr bool to_ulong_test()
{
    seqan3::dynamic_bitset t1{"0011000"};
    seqan3::dynamic_bitset t2{0b001100};

    return t1.to_ulong() == 24UL && t2.to_ulong() == 12UL;
}

TEST(dynamic_bitset, to_ulong)
{
    constexpr bool b = to_ulong_test();
    EXPECT_TRUE(b);
    EXPECT_TRUE(to_ulong_test());
    if constexpr (std::numeric_limits<unsigned long>::max() < std::numeric_limits<size_t>::max())
    {
        seqan3::dynamic_bitset t1{std::numeric_limits<unsigned long>::max() + 1};
        EXPECT_THROW(t1.to_ulong(), std::overflow_error);
    }
}

constexpr bool to_ullong_test()
{
    seqan3::dynamic_bitset t1{"0011000"};
    seqan3::dynamic_bitset t2{0b001100};

    return t1.to_ullong() == 24ULL && t2.to_ullong() == 12ULL;
}

TEST(dynamic_bitset, to_ullong)
{
    constexpr bool b = to_ullong_test();
    EXPECT_TRUE(b);
    EXPECT_TRUE(to_ullong_test());
    if constexpr (std::numeric_limits<unsigned long long>::max() < std::numeric_limits<size_t>::max())
    {
        seqan3::dynamic_bitset t1{std::numeric_limits<unsigned long long>::max() + 1};
        EXPECT_THROW(t1.to_ulong(), std::overflow_error);
    }
}

TEST(dynamic_bitset, output)
{
    seqan3::dynamic_bitset t1{"0011000"};
    std::ostringstream os;
    os << t1;
    EXPECT_EQ(os.str(), std::string{"0011000"});
}

TEST(dynamic_bitset, input)
{
    { // Until whitespace
        seqan3::dynamic_bitset t1{""};
        std::istringstream is{"0011 0001"};
        is >> t1;
        EXPECT_EQ(t1, seqan3::dynamic_bitset{"0011"});
    }

    { // Exceed capacity
        seqan3::dynamic_bitset<5> t1{"11111"};
        std::istringstream is{"00110001"};
        is >> t1;
        EXPECT_EQ(t1, seqan3::dynamic_bitset{"00110"});

        std::string remaining{};
        is >> remaining;
        EXPECT_EQ(remaining, std::string{"001"});
    }

    { // eof before capacity reached
        seqan3::dynamic_bitset t1{};
        std::istringstream is{"00110001"};
        is >> t1;
        EXPECT_EQ(t1, seqan3::dynamic_bitset{"00110001"});
    }
}

TEST(dynamic_bitset, debug_stream)
{
    std::ostringstream o;
    seqan3::debug_stream_type my_stream{o};

    seqan3::dynamic_bitset t1{0b1100'1110'1010'1111};

    my_stream << t1;
    o.flush();
    EXPECT_EQ(o.str(), "1100'1110'1010'1111");

    seqan3::dynamic_bitset const t2{0b1011'1010'1111'0000};

    my_stream << t2;
    o.flush();
    EXPECT_EQ(o.str(), "1100'1110'1010'11111011'1010'1111'0000");

    my_stream << seqan3::dynamic_bitset{0b0101'1110'0101'1001}; // The leftmost 0 will be stripped
    o.flush();
    EXPECT_EQ(o.str(), "1100'1110'1010'11111011'1010'1111'00001011'1100'1011'001");
}

TEST(dynamic_bitset, std_hash)
{
    seqan3::dynamic_bitset t1{"0011000"};
    seqan3::dynamic_bitset t2{0b001100};
    std::hash<seqan3::dynamic_bitset<58>> hasher{};

    EXPECT_EQ(hasher(t1), 24ULL);
    EXPECT_EQ(hasher(t2), 12ULL);
}

TEST(dynamic_bitset, serialisation)
{
    seqan3::dynamic_bitset t1{0b100101};
    seqan3::test::do_serialisation(t1);
}

// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <sstream>
#include <vector>

#include <gtest/gtest.h>

#include <seqan3/range/detail/random_access_iterator.hpp>

#include "../iterator_test_template.hpp"

// -----------------------------------------------------------------------------
// test templates
// -----------------------------------------------------------------------------

template <>
struct iterator_fixture<seqan3::detail::random_access_iterator<std::vector<int>>> : public ::testing::Test
{
    using iterator_tag = std::random_access_iterator_tag;
    static constexpr bool const_iterable = true;

    struct expose_iterator
    {
        using iterator_type = seqan3::detail::random_access_iterator<std::vector<int>>;
        using const_iterator_type = seqan3::detail::random_access_iterator<std::vector<int> const>;

        std::vector<int> rng{1, 2, 3, 4, 5, 6, 7, 8};

        iterator_type begin() { return iterator_type{rng}; }
        iterator_type end() { return iterator_type{rng, rng.size()}; }
        const_iterator_type begin() const { return const_iterator_type{rng}; }
        const_iterator_type end() const { return const_iterator_type{rng, rng.size()}; }
    };

    std::vector<int> expected_range{1, 2, 3, 4, 5, 6, 7, 8};
    expose_iterator test_range{};
};

using test_type = ::testing::Types<seqan3::detail::random_access_iterator<std::vector<int>>>;

INSTANTIATE_TYPED_TEST_SUITE_P(iterator_fixture, iterator_fixture, test_type, );

// -----------------------------------------------------------------------------
// individual tests
// -----------------------------------------------------------------------------

class random_access_iterator_test_fixture : public ::testing::Test
{
protected:
    std::vector<uint8_t> v_empty{};
    std::vector<uint8_t> const v_const_empty{};
    std::vector<uint8_t> v, v2, v3, v4, w, w2;
    std::vector<uint8_t> const v_const{'a', 't'};
    std::vector<uint8_t> const v2_const{'a', 'u'};
    std::vector<uint8_t> const v3_const{'a', 't', 'z'};
    std::vector<uint8_t> const v4_const{'a', 'u', 'v', 'w', 'x'};
    std::vector<uint8_t> const w_const{'c', 't'};
    std::vector<uint8_t> const w2_const{'b', 'v'};
    std::array<long int, 3> a;
    std::array<long int, 3> const a_const = {11, 22, 33};

    virtual void SetUp()
    {
        // code here will execute just before the test ensues
        v = {'a', 't'};
        v2 = {'a', 'u'};
        v3 = {'a', 't', 'z'};
        v4 = {'a', 'u', 'v', 'w', 'x'};
        w = {'c', 't'};
        w2 = {'b', 'v'};
        a = {11, 22, 33};
    }
};

// default constructor
TEST(random_access_iterator_test, default_constructor)
{
    [[maybe_unused]] seqan3::detail::random_access_iterator<std::vector<uint8_t>> it;
    [[maybe_unused]] seqan3::detail::random_access_iterator<std::vector<uint8_t> const> it2;
}

// constructor with empty container reference
TEST(random_access_iterator_test, constructor_ref)
{
    // non-const version
    std::vector<uint8_t> v_empty;
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it(v_empty);
    // const version
    std::vector<uint8_t> const v_const_empty;
    seqan3::detail::random_access_iterator<std::vector<uint8_t> const> it2(v_const_empty);
}

// test constructor call with non-empty container reference and subscript operator
TEST_F(random_access_iterator_test_fixture, constructor_ref2)
{
    // non-const version
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it(v);
    EXPECT_EQ('a', it[0]);
    EXPECT_EQ('t', it[1]);
    // const version
    seqan3::detail::random_access_iterator<std::vector<uint8_t> const> it2(v_const);
    EXPECT_EQ('a', it2[0]);
    EXPECT_EQ('t', it2[1]);
}

// test constructor with non-empty container reference and position offset  and subscript operator
TEST_F(random_access_iterator_test_fixture, constructor_ref3)
{
    // non-const version
    seqan3::detail::random_access_iterator<std::array<long int, 3>> it(a, 1);
    EXPECT_EQ(22, it[0]);
    EXPECT_EQ(33, it[1]);
    // const version
    seqan3::detail::random_access_iterator<std::array<long int, 3> const> it2(a_const, 1);
    EXPECT_EQ(22, it2[0]);
    EXPECT_EQ(33, it2[1]);
}

// copy constructor with empty container reference
TEST_F(random_access_iterator_test_fixture, cp_constructor1)
{
    // non-const container
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_base(v_empty);
    [[maybe_unused]] seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_derivate(it_base);
    // const container
    seqan3::detail::random_access_iterator<std::vector<uint8_t> const> it_base2(v_const_empty);
    [[maybe_unused]] seqan3::detail::random_access_iterator<std::vector<uint8_t> const> it_derivate2(it_base2);
}

// copy constructor with non-empty container reference
TEST_F(random_access_iterator_test_fixture, cp_constructor2)
{
    // non-const
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_base(v);
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_derivate(it_base);
    EXPECT_EQ('a', it_base[0]);
    EXPECT_EQ('a', it_derivate[0]);
    // const
    seqan3::detail::random_access_iterator<std::vector<uint8_t> const> it_base2(v_const);
    seqan3::detail::random_access_iterator<std::vector<uint8_t> const> it_derivate2(it_base2);
    EXPECT_EQ('a', it_base2[0]);
    EXPECT_EQ('a', it_derivate2[0]);

}

// construct const_iterator from iterator
TEST_F(random_access_iterator_test_fixture, constructor_const_from_nonconst)
{
    // non-const
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it(v);
    EXPECT_EQ('a', it[0]);
    // const
    seqan3::detail::random_access_iterator<std::vector<uint8_t> const> cit(it);
    EXPECT_EQ('a', cit[0]);
}

// test assignment construction with empty container reference
TEST_F(random_access_iterator_test_fixture, constructor_assign1)
{
    // non-const
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_base(v_empty);
    [[maybe_unused]] seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_derived = it_base;
    // const
    seqan3::detail::random_access_iterator<std::vector<uint8_t> const> it_base2(v_const_empty);
    [[maybe_unused]] seqan3::detail::random_access_iterator<std::vector<uint8_t> const> it_derived2 = it_base2;
}

// test assignment construction with non-empty container reference and subscript
TEST_F(random_access_iterator_test_fixture, constructor_assign2)
{
    // non-const
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_base(v);
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_derivate = it_base;
    EXPECT_EQ('t', it_base[1]);
    EXPECT_EQ('t', it_derivate[1]);
    // const
    seqan3::detail::random_access_iterator<std::vector<uint8_t> const> it_base2(v_const);
    seqan3::detail::random_access_iterator<std::vector<uint8_t> const> it_derivate2 = it_base2;
    EXPECT_EQ('t', it_base2[1]);
    EXPECT_EQ('t', it_derivate2[1]);
}

// test move constructor
TEST_F(random_access_iterator_test_fixture, constructor_move)
{
    // non-const
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it1(v);
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it2(std::move(it1));
    EXPECT_EQ('a', it2[0]);
    EXPECT_EQ('t', it2[1]);
    // const
    seqan3::detail::random_access_iterator<std::vector<uint8_t> const> it3(v_const);
    seqan3::detail::random_access_iterator<std::vector<uint8_t> const> it4(std::move(it3));
    EXPECT_EQ('a', it3[0]);
    EXPECT_EQ('t', it4[1]);
}

// test move assignment
TEST_F(random_access_iterator_test_fixture, move_assign)
{
    // non-const
    seqan3::detail::random_access_iterator<std::array<long int, 3>> it1, it2;
    it2 = std::move(it1);
    // const
    seqan3::detail::random_access_iterator<std::array<long int, 3> const> it3, it4;
    it4 = std::move(it3);
}

// explicit desctructor call
TEST_F(random_access_iterator_test_fixture, cp_destructor)
{
    // non-const
    using iterator_type = typename seqan3::detail::random_access_iterator<std::vector<uint8_t>>;
    iterator_type it(v_empty);
    iterator_type * it_ptr = &it;
    it_ptr->iterator_type::~iterator_type();
    // const
    using iterator_type2 = typename seqan3::detail::random_access_iterator<std::vector<uint8_t> const>;
    iterator_type2 it2(v_const_empty);
    iterator_type2 * it_ptr2 = &it2;
    it_ptr2->iterator_type2::~iterator_type2();

}

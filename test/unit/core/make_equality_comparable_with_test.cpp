// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2018, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/core/make_equality_comparable_with.hpp>
#include <seqan3/std/concepts>

template <typename type>
class make_equality_comparable_with_test : public ::testing::Test
{
public:

    auto get_test_value()
    {
        if constexpr (std::is_same_v<type, int32_t>)
        {
            return 10;
        }
        else if constexpr (std::is_same_v<type, float>)
        {
            return 3.2f;
        }
        else if constexpr (std::is_same_v<type, std::tuple<int32_t, int16_t>>)
        {
            return std::tuple{static_cast<int32_t>(10), static_cast<int16_t>(-1)};
        }
        else
        {
            FAIL(); //Did you added a new type?
        }
    }

    auto get_test_value_less()
    {
        if constexpr (std::is_same_v<type, int32_t>)
        {
            return 9;
        }
        else if constexpr (std::is_same_v<type, float>)
        {
            return 1.2f;
        }
        else if constexpr (std::is_same_v<type, std::tuple<int32_t, int16_t>>)
        {
            return std::tuple{static_cast<int32_t>(10), static_cast<int16_t>(-10)};
        }
        else
        {
            FAIL(); //Did you added a new type?
        }
    }

};

using test_types = testing::Types<int32_t, float, std::tuple<int32_t, int16_t>>;
TYPED_TEST_CASE(make_equality_comparable_with_test, test_types);

template <typename host_type>
struct host_wrapper
{
    host_type host;
};

template <typename host_type>
host_wrapper(host_type) -> host_wrapper<host_type>;

template <typename host_type>
class operator_tester :
    public seqan3::make_equality_comparable_with<operator_tester<host_type>, host_type>
{
public:

    operator_tester() = default;
    operator_tester(operator_tester const &) = default;
    operator_tester(operator_tester &&) = default;
    operator_tester & operator=(operator_tester const &) = default;
    operator_tester & operator=(operator_tester &&) = default;
    ~operator_tester() = default;

    explicit operator_tester(host_wrapper<host_type> const v) : host{v.host}
    {}

    constexpr operator host_type() const
    {
        return host;
    }

private:

    friend seqan3::make_equality_comparable_with<operator_tester<host_type>, host_type>;

    host_type const & compare_value() const
    {
        return host;
    }

    host_type host;
};

TYPED_TEST(make_equality_comparable_with_test, construct)
{
    using test_t = operator_tester<TypeParam>;
    EXPECT_TRUE(std::is_default_constructible_v<test_t>);
    EXPECT_TRUE(std::is_copy_constructible_v<test_t>);
    EXPECT_TRUE(std::is_move_constructible_v<test_t>);
    EXPECT_TRUE(std::is_copy_assignable_v<test_t>);
    EXPECT_TRUE(std::is_move_assignable_v<test_t>);
    EXPECT_TRUE(std::is_destructible_v<test_t>);
}

TYPED_TEST(make_equality_comparable_with_test, equality_concept)
{
    using test_t = operator_tester<TypeParam>;
    EXPECT_TRUE(std::EqualityComparable<test_t>);
    EXPECT_TRUE((std::EqualityComparableWith<test_t, TypeParam>));
}

TYPED_TEST(make_equality_comparable_with_test, cmp_eq)
{
    operator_tester<TypeParam> t1{host_wrapper{this->get_test_value()}};
    operator_tester<TypeParam> t2{host_wrapper{this->get_test_value()}};

    EXPECT_TRUE(t1 == t2);
    EXPECT_TRUE(t1 == this->get_test_value());
    EXPECT_TRUE(this->get_test_value() == t1);
}

TYPED_TEST(make_equality_comparable_with_test, cmp_ne)
{
    operator_tester<TypeParam> t1{host_wrapper{this->get_test_value()}};
    operator_tester<TypeParam> t2{host_wrapper{this->get_test_value()}};

    EXPECT_FALSE(t1 != t2);
    EXPECT_FALSE(t1 != this->get_test_value());
    EXPECT_FALSE(this->get_test_value() != t1);
}

TYPED_TEST(make_equality_comparable_with_test, cmp_lt)
{
    operator_tester<TypeParam> t1{host_wrapper{this->get_test_value()}};
    operator_tester<TypeParam> t2{host_wrapper{this->get_test_value_less()}};

    EXPECT_FALSE(t1 < t2);
    EXPECT_TRUE(t2 < t1);
    EXPECT_FALSE(t1 < this->get_test_value_less());
    EXPECT_TRUE(this->get_test_value_less() < t1);
}

TYPED_TEST(make_equality_comparable_with_test, cmp_le)
{
    operator_tester<TypeParam> t1{host_wrapper{this->get_test_value()}};
    operator_tester<TypeParam> t2{host_wrapper{this->get_test_value_less()}};

    EXPECT_FALSE(t1 <= t2);
    EXPECT_TRUE(t2 <= t1);
    EXPECT_TRUE(t1 <= t1);
    EXPECT_TRUE(t1 <= this->get_test_value());
    EXPECT_FALSE(t1 <= this->get_test_value_less());
    EXPECT_TRUE(this->get_test_value_less() <= t1);
    EXPECT_TRUE(this->get_test_value() <= t1);
}

TYPED_TEST(make_equality_comparable_with_test, cmp_gt)
{
    operator_tester<TypeParam> t1{host_wrapper{this->get_test_value()}};
    operator_tester<TypeParam> t2{host_wrapper{this->get_test_value_less()}};

    EXPECT_TRUE(t1 > t2);
    EXPECT_FALSE(t2 > t1);
    EXPECT_TRUE(t1 > this->get_test_value_less());
    EXPECT_FALSE(this->get_test_value_less() > t1);
}

TYPED_TEST(make_equality_comparable_with_test, cmp_ge)
{
    operator_tester<TypeParam> t1{host_wrapper{this->get_test_value()}};
    operator_tester<TypeParam> t2{host_wrapper{this->get_test_value_less()}};

    EXPECT_TRUE(t1 >= t2);
    EXPECT_FALSE(t2 >= t1);
    EXPECT_TRUE(t1 >= t1);
    EXPECT_TRUE(t1 >= this->get_test_value());
    EXPECT_TRUE(t1 >= this->get_test_value_less());
    EXPECT_FALSE(this->get_test_value_less() >= t1);
    EXPECT_TRUE(this->get_test_value() >= t1);
}

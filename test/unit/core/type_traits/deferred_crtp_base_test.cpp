// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <type_traits>

#include <seqan3/core/type_traits/deferred_crtp_base.hpp>

// Defines a crtp_base class with an additional value type.
template <typename derived_t, typename value_t = std::string>
class base1
{
public:

    value_t func1() const
    {
        return {"instance of base1"};
    }
};

// Defines a crtp_base class with an additional return type and a parameter_t.
template <typename derived_t, typename value_t = int, typename parameter_t = int>
class base2
{
public:

    value_t func2(parameter_t const p) const
    {
        return static_cast<value_t>(p);
    }
};

// The derived class that inherits from a variadic crtp pattern.
template <typename ...bases_t>
class derived : public seqan3::detail::invoke_deferred_crtp_base<bases_t, derived<bases_t...>>...
{};

TEST(deferred_crtp_base, one_base_not_augmented)
{
    using deferred_base1 = seqan3::detail::deferred_crtp_base<base1>;

    derived<deferred_base1> d{};

    EXPECT_TRUE((std::is_same_v<decltype(d.func1()), std::string>));
}

TEST(deferred_crtp_base, multiple_base_not_augmented)
{
    using deferred_base1 = seqan3::detail::deferred_crtp_base<base1>;
    using deferred_base2 = seqan3::detail::deferred_crtp_base<base2>;

    derived<deferred_base1, deferred_base2> d{};

    EXPECT_TRUE((std::is_same_v<decltype(d.func1()), std::string>));
    EXPECT_TRUE((std::is_same_v<decltype(d.func2(10)), int>));
}

TEST(deferred_crtp_base, one_base_augmented)
{
    using deferred_base1 = seqan3::detail::deferred_crtp_base<base1, std::vector<char>>;

    derived<deferred_base1> d{};

    EXPECT_TRUE((std::is_same_v<decltype(d.func1()), std::vector<char>>));
}

TEST(deferred_crtp_base, multiple_base_augmented)
{
    using deferred_base1 = seqan3::detail::deferred_crtp_base<base1, std::vector<char>>;
    using deferred_base2 = seqan3::detail::deferred_crtp_base<base2, float, int8_t>;

    derived<deferred_base1, deferred_base2> d{};

    EXPECT_TRUE((std::is_same_v<decltype(d.func1()), std::vector<char>>));
    EXPECT_TRUE((std::is_same_v<decltype(d.func2(10)), float>));
}

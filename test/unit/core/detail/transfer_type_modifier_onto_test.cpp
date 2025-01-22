// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/core/detail/transfer_type_modifier_onto.hpp>
#include <seqan3/test/expect_same_type.hpp>

TEST(transfer_type_modifier_onto, type)
{
    // target type has no modifier
    EXPECT_SAME_TYPE(double, (typename seqan3::detail::transfer_type_modifier_onto<int, double>::type));
    EXPECT_SAME_TYPE(double &, (typename seqan3::detail::transfer_type_modifier_onto<int &, double>::type));
    EXPECT_SAME_TYPE(double &&, (typename seqan3::detail::transfer_type_modifier_onto<int &&, double>::type));
    EXPECT_SAME_TYPE(double const, (typename seqan3::detail::transfer_type_modifier_onto<int const, double>::type));
    EXPECT_SAME_TYPE(double const &, (typename seqan3::detail::transfer_type_modifier_onto<int const &, double>::type));
    EXPECT_SAME_TYPE(double const &&,
                     (typename seqan3::detail::transfer_type_modifier_onto<int const &&, double>::type));

    // target type has lvalue reference modifier
    EXPECT_SAME_TYPE(double &, (typename seqan3::detail::transfer_type_modifier_onto<int, double &>::type));
    EXPECT_SAME_TYPE(double & /*&*/, (typename seqan3::detail::transfer_type_modifier_onto<int &, double &>::type));
    EXPECT_SAME_TYPE(double & /*&&*/, (typename seqan3::detail::transfer_type_modifier_onto<int &&, double &>::type));
    EXPECT_SAME_TYPE(double const &, (typename seqan3::detail::transfer_type_modifier_onto<int const, double &>::type));
    EXPECT_SAME_TYPE(double const & /*&*/,
                     (typename seqan3::detail::transfer_type_modifier_onto<int const &, double &>::type));
    EXPECT_SAME_TYPE(double const & /*&&*/,
                     (typename seqan3::detail::transfer_type_modifier_onto<int const &&, double &>::type));

    // target type has rvalue reference modifier
    EXPECT_SAME_TYPE(double &&, (typename seqan3::detail::transfer_type_modifier_onto<int, double &&>::type));
    EXPECT_SAME_TYPE(double /*&&*/ &, (typename seqan3::detail::transfer_type_modifier_onto<int &, double &&>::type));
    EXPECT_SAME_TYPE(double && /*&&*/, (typename seqan3::detail::transfer_type_modifier_onto<int &&, double &&>::type));
    EXPECT_SAME_TYPE(double const &&,
                     (typename seqan3::detail::transfer_type_modifier_onto<int const, double &&>::type));
    EXPECT_SAME_TYPE(double const /*&&*/ &,
                     (typename seqan3::detail::transfer_type_modifier_onto<int const &, double &&>::type));
    EXPECT_SAME_TYPE(double const /*&&*/ &&,
                     (typename seqan3::detail::transfer_type_modifier_onto<int const &&, double &&>::type));

    // target type has const modifier
    EXPECT_SAME_TYPE(double const, (typename seqan3::detail::transfer_type_modifier_onto<int, double const>::type));
    EXPECT_SAME_TYPE(double const &, (typename seqan3::detail::transfer_type_modifier_onto<int &, double const>::type));
    EXPECT_SAME_TYPE(double const &&,
                     (typename seqan3::detail::transfer_type_modifier_onto<int &&, double const>::type));
    EXPECT_SAME_TYPE(double /*const*/ const,
                     (typename seqan3::detail::transfer_type_modifier_onto<int const, double const>::type));
    EXPECT_SAME_TYPE(double /*const*/ const &,
                     (typename seqan3::detail::transfer_type_modifier_onto<int const &, double const>::type));
    EXPECT_SAME_TYPE(double /*const*/ const &&,
                     (typename seqan3::detail::transfer_type_modifier_onto<int const &&, double const>::type));

    // target type has const lvalue reference modifier
    EXPECT_SAME_TYPE(double const &, (typename seqan3::detail::transfer_type_modifier_onto<int, double const &>::type));
    EXPECT_SAME_TYPE(double const & /*&*/,
                     (typename seqan3::detail::transfer_type_modifier_onto<int &, double const &>::type));
    EXPECT_SAME_TYPE(double const & /*&&*/,
                     (typename seqan3::detail::transfer_type_modifier_onto<int &&, double const &>::type));
    EXPECT_SAME_TYPE(double /*const*/ const &,
                     (typename seqan3::detail::transfer_type_modifier_onto<int const, double const &>::type));
    EXPECT_SAME_TYPE(double /*const*/ const & /*&*/,
                     (typename seqan3::detail::transfer_type_modifier_onto<int const &, double const &>::type));
    EXPECT_SAME_TYPE(double /*const*/ const & /*&&*/,
                     (typename seqan3::detail::transfer_type_modifier_onto<int const &&, double const &>::type));

    // target type has const rvalue reference modifier
    EXPECT_SAME_TYPE(double const &&,
                     (typename seqan3::detail::transfer_type_modifier_onto<int, double const &&>::type));
    EXPECT_SAME_TYPE(double const /*&&*/ &,
                     (typename seqan3::detail::transfer_type_modifier_onto<int &, double const &&>::type));
    EXPECT_SAME_TYPE(double const && /*&&*/,
                     (typename seqan3::detail::transfer_type_modifier_onto<int &&, double const &&>::type));
    EXPECT_SAME_TYPE(double /*const*/ const &&,
                     (typename seqan3::detail::transfer_type_modifier_onto<int const, double const &&>::type));
    EXPECT_SAME_TYPE(double /*const*/ const /*&&*/ &,
                     (typename seqan3::detail::transfer_type_modifier_onto<int const &, double const &&>::type));
    EXPECT_SAME_TYPE(double /*const*/ const /*&&*/ &&,
                     (typename seqan3::detail::transfer_type_modifier_onto<int const &&, double const &&>::type));
}

TEST(transfer_type_modifier_onto, type_t_helper)
{
    // target type has no modifier
    EXPECT_SAME_TYPE(double, (seqan3::detail::transfer_type_modifier_onto_t<int, double>));
    EXPECT_SAME_TYPE(double &, (seqan3::detail::transfer_type_modifier_onto_t<int &, double>));
    EXPECT_SAME_TYPE(double &&, (seqan3::detail::transfer_type_modifier_onto_t<int &&, double>));
    EXPECT_SAME_TYPE(double const, (seqan3::detail::transfer_type_modifier_onto_t<int const, double>));
    EXPECT_SAME_TYPE(double const &, (seqan3::detail::transfer_type_modifier_onto_t<int const &, double>));
    EXPECT_SAME_TYPE(double const &&, (seqan3::detail::transfer_type_modifier_onto_t<int const &&, double>));

    // target type has lvalue reference modifier
    EXPECT_SAME_TYPE(double &, (seqan3::detail::transfer_type_modifier_onto_t<int, double &>));
    EXPECT_SAME_TYPE(double & /*&*/, (seqan3::detail::transfer_type_modifier_onto_t<int &, double &>));
    EXPECT_SAME_TYPE(double & /*&&*/, (seqan3::detail::transfer_type_modifier_onto_t<int &&, double &>));
    EXPECT_SAME_TYPE(double const &, (seqan3::detail::transfer_type_modifier_onto_t<int const, double &>));
    EXPECT_SAME_TYPE(double const & /*&*/, (seqan3::detail::transfer_type_modifier_onto_t<int const &, double &>));
    EXPECT_SAME_TYPE(double const & /*&&*/, (seqan3::detail::transfer_type_modifier_onto_t<int const &&, double &>));

    // target type has rvalue reference modifier
    EXPECT_SAME_TYPE(double &&, (seqan3::detail::transfer_type_modifier_onto_t<int, double &&>));
    EXPECT_SAME_TYPE(double /*&&*/ &, (seqan3::detail::transfer_type_modifier_onto_t<int &, double &&>));
    EXPECT_SAME_TYPE(double && /*&&*/, (seqan3::detail::transfer_type_modifier_onto_t<int &&, double &&>));
    EXPECT_SAME_TYPE(double const &&, (seqan3::detail::transfer_type_modifier_onto_t<int const, double &&>));
    EXPECT_SAME_TYPE(double const /*&&*/ &, (seqan3::detail::transfer_type_modifier_onto_t<int const &, double &&>));
    EXPECT_SAME_TYPE(double const /*&&*/ &&, (seqan3::detail::transfer_type_modifier_onto_t<int const &&, double &&>));

    // target type has const modifier
    EXPECT_SAME_TYPE(double const, (seqan3::detail::transfer_type_modifier_onto_t<int, double const>));
    EXPECT_SAME_TYPE(double const &, (seqan3::detail::transfer_type_modifier_onto_t<int &, double const>));
    EXPECT_SAME_TYPE(double const &&, (seqan3::detail::transfer_type_modifier_onto_t<int &&, double const>));
    EXPECT_SAME_TYPE(double /*const*/ const, (seqan3::detail::transfer_type_modifier_onto_t<int const, double const>));
    EXPECT_SAME_TYPE(double /*const*/ const &,
                     (seqan3::detail::transfer_type_modifier_onto_t<int const &, double const>));
    EXPECT_SAME_TYPE(double /*const*/ const &&,
                     (seqan3::detail::transfer_type_modifier_onto_t<int const &&, double const>));

    // target type has const lvalue reference modifier
    EXPECT_SAME_TYPE(double const &, (seqan3::detail::transfer_type_modifier_onto_t<int, double const &>));
    EXPECT_SAME_TYPE(double const & /*&*/, (seqan3::detail::transfer_type_modifier_onto_t<int &, double const &>));
    EXPECT_SAME_TYPE(double const & /*&&*/, (seqan3::detail::transfer_type_modifier_onto_t<int &&, double const &>));
    EXPECT_SAME_TYPE(double /*const*/ const &,
                     (seqan3::detail::transfer_type_modifier_onto_t<int const, double const &>));
    EXPECT_SAME_TYPE(double /*const*/ const & /*&*/,
                     (seqan3::detail::transfer_type_modifier_onto_t<int const &, double const &>));
    EXPECT_SAME_TYPE(double /*const*/ const & /*&&*/,
                     (seqan3::detail::transfer_type_modifier_onto_t<int const &&, double const &>));

    // target type has const rvalue reference modifier
    EXPECT_SAME_TYPE(double const &&, (seqan3::detail::transfer_type_modifier_onto_t<int, double const &&>));
    EXPECT_SAME_TYPE(double const /*&&*/ &, (seqan3::detail::transfer_type_modifier_onto_t<int &, double const &&>));
    EXPECT_SAME_TYPE(double const && /*&&*/, (seqan3::detail::transfer_type_modifier_onto_t<int &&, double const &&>));
    EXPECT_SAME_TYPE(double /*const*/ const &&,
                     (seqan3::detail::transfer_type_modifier_onto_t<int const, double const &&>));
    EXPECT_SAME_TYPE(double /*const*/ const /*&&*/ &,
                     (seqan3::detail::transfer_type_modifier_onto_t<int const &, double const &&>));
    EXPECT_SAME_TYPE(double /*const*/ const /*&&*/ &&,
                     (seqan3::detail::transfer_type_modifier_onto_t<int const &&, double const &&>));
}

// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <concepts>
#include <list>
#include <vector>

#include <seqan3/utility/type_traits/lazy_conditional.hpp>

TEST(lazy, instantiate)
{
    EXPECT_TRUE((std::is_same_v<seqan3::detail::instantiate_t<std::vector<int>>, std::vector<int>>));

    EXPECT_TRUE(
        (std::is_same_v<seqan3::detail::instantiate_t<seqan3::detail::lazy<std::vector, int>>, std::vector<int>>));
}

template <typename t>
    requires std::integral<t>
using integral_identity_t = t;

TEST(lazy, lazy_conditional)
{
    // regular conditional behaviour
    EXPECT_TRUE(
        (std::is_same_v<seqan3::detail::lazy_conditional_t<true, std::true_type, std::false_type>, std::true_type>));
    EXPECT_TRUE(
        (std::is_same_v<seqan3::detail::lazy_conditional_t<false, std::true_type, std::false_type>, std::false_type>));

    // lazy behaviour, safe
    EXPECT_TRUE((std::is_same_v<seqan3::detail::lazy_conditional_t<true,
                                                                   seqan3::detail::lazy<std::vector, int>,
                                                                   seqan3::detail::lazy<std::list, int>>,
                                std::vector<int>>));
    EXPECT_TRUE((std::is_same_v<seqan3::detail::lazy_conditional_t<false,
                                                                   seqan3::detail::lazy<std::vector, int>,
                                                                   seqan3::detail::lazy<std::list, int>>,
                                std::list<int>>));

    // lazy behaviour, important
    EXPECT_TRUE(
        (std::is_same_v<seqan3::detail::lazy_conditional_t<true, seqan3::detail::lazy<integral_identity_t, int>, void>,
                        int>));
    EXPECT_TRUE(
        (std::is_same_v<seqan3::detail::lazy_conditional_t<false, void, seqan3::detail::lazy<integral_identity_t, int>>,
                        int>));
}

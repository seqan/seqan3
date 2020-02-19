// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <list>
#include <vector>

#include <seqan3/core/type_traits/lazy.hpp>
#include <seqan3/std/concepts>

TEST(lazy, instantiate)
{
    EXPECT_TRUE((std::is_same_v<seqan3::detail::instantiate_t<std::vector<int>>,
                                std::vector<int>>));

    EXPECT_TRUE((std::is_same_v<seqan3::detail::instantiate_t<seqan3::detail::lazy<std::vector, int>>,
                                std::vector<int>>));
}

template <typename t>
    requires std::integral<t>
using integral_identity_t = t;

TEST(lazy, lazy_conditional)
{
    // regular conditional behaviour
    EXPECT_TRUE((std::is_same_v<seqan3::detail::lazy_conditional_t<true,
                                                                   std::true_type,
                                                                   std::false_type>,
                                std::true_type>));
    EXPECT_TRUE((std::is_same_v<seqan3::detail::lazy_conditional_t<false,
                                                                   std::true_type,
                                                                   std::false_type>,
                                std::false_type>));

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
    EXPECT_TRUE((std::is_same_v<seqan3::detail::lazy_conditional_t<true,
                                                                   seqan3::detail::lazy<integral_identity_t, int>,
                                                                   void>,
                                int>));
    EXPECT_TRUE((std::is_same_v<seqan3::detail::lazy_conditional_t<false,
                                                                   void,
                                                                   seqan3::detail::lazy<integral_identity_t, int>>,
                                int>));
}

// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <type_traits>

#include <seqan3/utility/type_list/type_list.hpp>

TEST(type_list, basic)
{
    using t = seqan3::type_list<int, char, double>;

    EXPECT_TRUE((std::is_same_v<typename t::type, t>));
    EXPECT_EQ(t::size(), 3u);
}

TEST(type_list, concept)
{
    using t = seqan3::type_list<int, char, double>;

    EXPECT_TRUE((seqan3::detail::template_specialisation_of<t, seqan3::type_list>));
    EXPECT_FALSE((seqan3::detail::template_specialisation_of<int, seqan3::type_list>));
}

// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

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

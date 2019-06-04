// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <list>
#include <vector>

#include <gtest/gtest.h>

#include <meta/meta.hpp>

#include <range/v3/view/iota.hpp>
#include <range/v3/view/take_exactly.hpp>

#include <seqan3/core/type_traits/basic.hpp>
#include <seqan3/range/detail/random_access_iterator.hpp>

using namespace seqan3;

TEST(type_trait, remove_cvref_t)
{
    EXPECT_TRUE((std::is_same_v<int, remove_cvref_t<int>>));
    EXPECT_TRUE((std::is_same_v<int, remove_cvref_t<int const>>));
    EXPECT_TRUE((std::is_same_v<int, remove_cvref_t<int volatile>>));
    EXPECT_TRUE((std::is_same_v<int, remove_cvref_t<int &>>));
    EXPECT_TRUE((std::is_same_v<int, remove_cvref_t<int &&>>));
    EXPECT_TRUE((std::is_same_v<int, remove_cvref_t<int const &>>));
    EXPECT_TRUE((std::is_same_v<int, remove_cvref_t<int const &&>>));
    EXPECT_TRUE((std::is_same_v<int, remove_cvref_t<int const volatile &&>>));
    // don't decay pointers and arrays:
    EXPECT_FALSE((std::is_same_v<int, remove_cvref_t<int*>>));      // type stays same
    EXPECT_FALSE((std::is_same_v<int, remove_cvref_t<int[3]>>));    // type stays same
    EXPECT_FALSE((std::is_same_v<int*, remove_cvref_t<int[3]>>));   // type stays same
    // the last example would be true for std::decay_t
}

// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <iostream>
#include <memory>

#include <seqan3/core/range/detail/adaptor_for_view_without_args.hpp>
#include <seqan3/utility/type_traits/basic.hpp>

template <typename t>
struct dummy_view
{};

TEST(adaptor_combination, constexpr_combine)
{
    constexpr auto adaptor1 = seqan3::detail::adaptor_for_view_without_args<dummy_view>{};
    constexpr auto adaptor2 = seqan3::detail::adaptor_for_view_without_args<dummy_view>{};

    EXPECT_TRUE((SEQAN3_IS_CONSTEXPR(adaptor1 | adaptor2)));
}

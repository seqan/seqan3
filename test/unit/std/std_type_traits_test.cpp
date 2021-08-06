// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/std/type_traits>

#include <seqan3/test/expect_same_type.hpp>

template <typename t>
SEQAN3_CONCEPT has_type = requires()
{
    typename t::type;
};

TEST(common_reference, zero_template_arguments)
{
    using common_reference = std::common_reference<>;
    EXPECT_FALSE(has_type<common_reference>);
}

TEST(common_reference, one_template_argument)
{
    EXPECT_SAME_TYPE(std::common_reference<int>::type,
                     int);
    EXPECT_SAME_TYPE(std::common_reference_t<int>,
                     int);
}

TEST(common_reference, two_template_argument_both_lvalue_reference_type)
{
    EXPECT_SAME_TYPE((std::common_reference<int &, long &>::type),
                     long);
    EXPECT_SAME_TYPE((std::common_reference_t<int &, long &>),
                     long);
}

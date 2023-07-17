// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/test/concept_helper_classes.hpp>
#include <seqan3/utility/concept.hpp>

// Requires that two types can be implicitly converted. Same as https://en.cppreference.com/w/cpp/types/is_convertible
TEST(convertability_concepts, implicitly_convertible_to)
{
    EXPECT_TRUE((seqan3::implicitly_convertible_to<type_a, type_a>));
    EXPECT_FALSE((seqan3::implicitly_convertible_to<type_a, type_b>));
    EXPECT_FALSE((seqan3::implicitly_convertible_to<type_a, type_c>));
    EXPECT_FALSE((seqan3::implicitly_convertible_to<type_a, type_d>));

    EXPECT_TRUE((seqan3::implicitly_convertible_to<type_b, type_a>)); // type_a is base of type_b
    EXPECT_TRUE((seqan3::implicitly_convertible_to<type_b, type_b>));
    EXPECT_TRUE((seqan3::implicitly_convertible_to<type_b, type_c>)); // custom constructor
    EXPECT_FALSE((seqan3::implicitly_convertible_to<type_b, type_d>));

    EXPECT_FALSE((seqan3::implicitly_convertible_to<type_c, type_a>));
    EXPECT_FALSE((seqan3::implicitly_convertible_to<type_c, type_b>));
    EXPECT_TRUE((seqan3::implicitly_convertible_to<type_c, type_c>));
    EXPECT_FALSE((seqan3::implicitly_convertible_to<type_c, type_d>));

    EXPECT_TRUE((seqan3::implicitly_convertible_to<type_d, type_a>));  // type_a is base of type_b is base of type_d
    EXPECT_TRUE((seqan3::implicitly_convertible_to<type_d, type_b>));  // type_b is base of type_d
    EXPECT_TRUE((seqan3::implicitly_convertible_to<type_d, type_c>));  // type_d -> type_a, then custom constructor
    EXPECT_FALSE((seqan3::implicitly_convertible_to<type_d, type_d>)); // unconstructible
}

// Requires that two types can be implicitly converted: static_cast<to>(from).
TEST(convertability_concepts, explicitly_convertible_to)
{
    EXPECT_TRUE((seqan3::explicitly_convertible_to<type_a, type_a>));
    EXPECT_FALSE((seqan3::explicitly_convertible_to<type_a, type_b>));
    EXPECT_TRUE((seqan3::explicitly_convertible_to<type_a, type_c>)); // custom constructor
    EXPECT_FALSE((seqan3::explicitly_convertible_to<type_a, type_d>));

    EXPECT_TRUE((seqan3::explicitly_convertible_to<type_b, type_a>)); // type_a is base of type_b
    EXPECT_TRUE((seqan3::explicitly_convertible_to<type_b, type_b>));
    EXPECT_TRUE((seqan3::explicitly_convertible_to<type_b, type_c>)); // custom constructor
    EXPECT_FALSE((seqan3::explicitly_convertible_to<type_b, type_d>));

    EXPECT_FALSE((seqan3::explicitly_convertible_to<type_c, type_a>));
    EXPECT_FALSE((seqan3::explicitly_convertible_to<type_c, type_b>));
    EXPECT_TRUE((seqan3::explicitly_convertible_to<type_c, type_c>));
    EXPECT_FALSE((seqan3::explicitly_convertible_to<type_c, type_d>));

    EXPECT_TRUE((seqan3::explicitly_convertible_to<type_d, type_a>));  // type_a is base of type_b is base of type_d
    EXPECT_TRUE((seqan3::explicitly_convertible_to<type_d, type_b>));  // type_b is base of type_d
    EXPECT_TRUE((seqan3::explicitly_convertible_to<type_d, type_c>));  // type_d -> type_a, then custom constructor
    EXPECT_FALSE((seqan3::explicitly_convertible_to<type_d, type_d>)); // unconstructible
}

// Requires that two types are both implicitly and explicitly convertible.
TEST(convertability_concepts, convertible_to)
{
    EXPECT_TRUE((std::convertible_to<type_a, type_a>));
    EXPECT_FALSE((std::convertible_to<type_a, type_b>));
    EXPECT_FALSE((std::convertible_to<type_a, type_c>));
    EXPECT_FALSE((std::convertible_to<type_a, type_d>));

    EXPECT_TRUE((std::convertible_to<type_b, type_a>));
    EXPECT_TRUE((std::convertible_to<type_b, type_b>));
    EXPECT_TRUE((std::convertible_to<type_b, type_c>));
    EXPECT_FALSE((std::convertible_to<type_b, type_d>));

    EXPECT_FALSE((std::convertible_to<type_c, type_a>));
    EXPECT_FALSE((std::convertible_to<type_c, type_b>));
    EXPECT_TRUE((std::convertible_to<type_c, type_c>));
    EXPECT_FALSE((std::convertible_to<type_c, type_d>));

    EXPECT_TRUE((std::convertible_to<type_d, type_a>));
    EXPECT_TRUE((std::convertible_to<type_d, type_b>));
    EXPECT_TRUE((std::convertible_to<type_d, type_c>));
    EXPECT_FALSE((std::convertible_to<type_d, type_d>));
}

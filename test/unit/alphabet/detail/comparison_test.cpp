// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/alphabet/detail/concept.hpp>
#include <seqan3/test/concept_helper_classes.hpp>

// Requires that two types can be compared for equality with each other (in either order) using both `==` and `!=`.
TEST(comparison_concepts, weakly_equality_comparable_with)
{
    EXPECT_FALSE((seqan3::detail::weakly_equality_comparable_with<type_a, type_a>));
    EXPECT_TRUE((seqan3::detail::weakly_equality_comparable_with<type_a, type_b>));
    EXPECT_FALSE((seqan3::detail::weakly_equality_comparable_with<type_a, type_c>));
    EXPECT_TRUE((seqan3::detail::weakly_equality_comparable_with<type_a, type_d>)); // via inheritance

    EXPECT_TRUE((seqan3::detail::weakly_equality_comparable_with<type_b, type_b>));
    EXPECT_TRUE((seqan3::detail::weakly_equality_comparable_with<type_b, type_c>));
    EXPECT_TRUE((seqan3::detail::weakly_equality_comparable_with<type_b, type_d>));

    EXPECT_TRUE((seqan3::detail::weakly_equality_comparable_with<type_c, type_c>));
    EXPECT_TRUE((seqan3::detail::weakly_equality_comparable_with<type_c, type_d>));

    EXPECT_TRUE((seqan3::detail::weakly_equality_comparable_with<type_d, type_d>));
}

// Requires that two types can be compared for equality with each other (in either order) using both `==` and `!=`.
// Additionally, both types need to be equality_comparible (with itself).
TEST(comparison_concepts, equality_comparable_with)
{
    EXPECT_FALSE((std::equality_comparable_with<type_a, type_a>));
    EXPECT_FALSE((std::equality_comparable_with<type_a, type_b>)); // type_a is not equality_comparable
    EXPECT_FALSE((std::equality_comparable_with<type_a, type_c>));
    EXPECT_FALSE((std::equality_comparable_with<type_a, type_d>)); // type_a is not equality_comparable

    EXPECT_TRUE((std::equality_comparable_with<type_b, type_b>));
    EXPECT_TRUE((std::equality_comparable_with<type_b, type_c>));
    EXPECT_TRUE((std::equality_comparable_with<type_b, type_d>));

    EXPECT_TRUE((std::equality_comparable_with<type_c, type_c>));
    EXPECT_TRUE((std::equality_comparable_with<type_c, type_d>));

    EXPECT_TRUE((std::equality_comparable_with<type_d, type_d>));
}

// Requires that two types can be compared in a partial order with each other (in either order) using <, >, <=, and >=.
// Also called `__PartiallyOrderedWith` in the STL (exposition only concept).
TEST(comparison_concepts, weakly_ordered_with)
{
    EXPECT_TRUE((seqan3::detail::weakly_ordered_with<type_a, type_a>));
    EXPECT_TRUE((seqan3::detail::weakly_ordered_with<type_a, type_b>));
    EXPECT_FALSE((seqan3::detail::weakly_ordered_with<type_a, type_c>));
    EXPECT_TRUE((seqan3::detail::weakly_ordered_with<type_a, type_d>)); // via inheritance

    EXPECT_TRUE((seqan3::detail::weakly_ordered_with<type_b, type_b>));
    EXPECT_FALSE((seqan3::detail::weakly_ordered_with<type_b, type_c>));
    EXPECT_TRUE((seqan3::detail::weakly_ordered_with<type_b, type_d>));

    EXPECT_FALSE((seqan3::detail::weakly_ordered_with<type_c, type_c>));
    EXPECT_FALSE((seqan3::detail::weakly_ordered_with<type_c, type_d>));

    EXPECT_TRUE((seqan3::detail::weakly_ordered_with<type_d, type_d>));
}

// Requires that two types can be compared order with each other (in either order) using ==, !=, <, >, <=, and >=.
TEST(comparison_concepts, totally_ordered_with)
{
    EXPECT_FALSE((std::totally_ordered_with<type_a, type_a>));
    EXPECT_FALSE((std::totally_ordered_with<type_a, type_b>)); // type_a is not totally_ordered
    EXPECT_FALSE((std::totally_ordered_with<type_a, type_c>));
    EXPECT_FALSE((std::totally_ordered_with<type_a, type_d>)); // type_a is not totally_ordered

    EXPECT_TRUE((std::totally_ordered_with<type_b, type_b>));
    EXPECT_FALSE((std::totally_ordered_with<type_b, type_c>));
    EXPECT_TRUE((std::totally_ordered_with<type_b, type_d>));

    EXPECT_FALSE((std::totally_ordered_with<type_c, type_c>));
    EXPECT_FALSE((std::totally_ordered_with<type_c, type_d>));

    EXPECT_TRUE((std::totally_ordered_with<type_d, type_d>));
}

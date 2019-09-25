// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <type_traits>

#include <range/v3/view/iota.hpp>

#include <seqan3/alignment/matrix/alignment_coordinate.hpp>
#include <seqan3/test/pretty_printing.hpp>

using namespace seqan3;

TEST(advanceable_alignment_coordinate, column_index_type)
{
    detail::column_index_type ci{1u};
    EXPECT_EQ(ci.get(), 1u);
    EXPECT_TRUE((std::same_as<std::remove_reference_t<decltype(ci.get())>, size_t>));

    detail::column_index_type ci2{1};
    EXPECT_EQ(ci2.get(), 1);
    EXPECT_TRUE((std::same_as<std::remove_reference_t<decltype(ci2.get())>, std::ptrdiff_t>));
}

TEST(advanceable_alignment_coordinate, row_index_type)
{
    detail::row_index_type ri{1u};
    EXPECT_EQ(ri.get(), 1u);
    EXPECT_TRUE((std::same_as<std::remove_reference_t<decltype(ri.get())>, size_t>));

    detail::row_index_type ri2{1};
    EXPECT_EQ(ri2.get(), 1);
    EXPECT_TRUE((std::same_as<std::remove_reference_t<decltype(ri2.get())>, std::ptrdiff_t>));
}

TEST(advanceable_alignment_coordinate, construction)
{
    EXPECT_TRUE(std::is_default_constructible<detail::advanceable_alignment_coordinate<>>::value);
    EXPECT_TRUE(std::is_copy_constructible<detail::advanceable_alignment_coordinate<>>::value);
    EXPECT_TRUE(std::is_copy_assignable<detail::advanceable_alignment_coordinate<>>::value);
    EXPECT_TRUE(std::is_move_constructible<detail::advanceable_alignment_coordinate<>>::value);
    EXPECT_TRUE(std::is_move_assignable<detail::advanceable_alignment_coordinate<>>::value);
    EXPECT_TRUE(std::is_destructible<detail::advanceable_alignment_coordinate<>>::value);
}

TEST(advanceable_alignment_coordinate, construction_with_different_state)
{
    detail::advanceable_alignment_coordinate<detail::advanceable_alignment_coordinate_state::row>
        ro{detail::column_index_type{2u}, detail::row_index_type{3u}};

    detail::advanceable_alignment_coordinate<detail::advanceable_alignment_coordinate_state::none> no{ro};

    EXPECT_EQ(no.first, 2u);
    EXPECT_EQ(no.second, 3u);
}

TEST(advanceable_alignment_coordinate, type_deduction)
{
    detail::advanceable_alignment_coordinate def_co{};
    EXPECT_TRUE((std::is_same_v<decltype(def_co),
                    detail::advanceable_alignment_coordinate<detail::advanceable_alignment_coordinate_state::none>>));

    detail::advanceable_alignment_coordinate co{detail::column_index_type{2u}, detail::row_index_type{3u}};
    EXPECT_TRUE((std::is_same_v<decltype(co),
                    detail::advanceable_alignment_coordinate<detail::advanceable_alignment_coordinate_state::none>>));
}

TEST(advanceable_alignment_coordinate, access)
{
    detail::advanceable_alignment_coordinate def_co{};
    EXPECT_EQ(def_co.first, 0u);
    EXPECT_EQ(def_co.second, 0u);

    detail::advanceable_alignment_coordinate co{detail::column_index_type{2u}, detail::row_index_type{3u}};
    EXPECT_EQ(co.first, 2u);
    EXPECT_EQ(co.second, 3u);
}

TEST(advanceable_alignment_coordinate, weakly_equality_comparable_concept)
{
    using not_incrementable =
        detail::advanceable_alignment_coordinate<detail::advanceable_alignment_coordinate_state::none>;
    using row_incrementable =
        detail::advanceable_alignment_coordinate<detail::advanceable_alignment_coordinate_state::row>;
    using column_incrementable =
        detail::advanceable_alignment_coordinate<detail::advanceable_alignment_coordinate_state::column>;

    EXPECT_TRUE(std::equality_comparable<not_incrementable>);
    EXPECT_TRUE(std::equality_comparable<row_incrementable>);
    EXPECT_TRUE(std::equality_comparable<column_incrementable>);
}

TEST(advanceable_alignment_coordinate, equality)
{
    using test_type =
        detail::advanceable_alignment_coordinate<detail::advanceable_alignment_coordinate_state::none>;

    test_type t1{detail::column_index_type{10u}, detail::row_index_type{5u}};
    test_type t2{detail::column_index_type{5u}, detail::row_index_type{5u}};
    test_type t3{detail::column_index_type{10u}, detail::row_index_type{10u}};

    EXPECT_TRUE(t1 == t1);
    EXPECT_FALSE(t2 == t1);
    EXPECT_FALSE(t1 == t3);
    EXPECT_FALSE(t2 == t3);
}

TEST(advanceable_alignment_coordinate, inequality)
{
    using test_type =
        detail::advanceable_alignment_coordinate<detail::advanceable_alignment_coordinate_state::none>;

    test_type t1{detail::column_index_type{10u}, detail::row_index_type{5u}};
    test_type t2{detail::column_index_type{5u}, detail::row_index_type{5u}};
    test_type t3{detail::column_index_type{10u}, detail::row_index_type{10u}};

    EXPECT_FALSE(t1 != t1);
    EXPECT_TRUE(t2 != t1);
    EXPECT_TRUE(t1 != t3);
    EXPECT_TRUE(t2 != t3);
}

TEST(advanceable_alignment_coordinate, incremental_concept)
{
    using not_incrementable =
        detail::advanceable_alignment_coordinate<detail::advanceable_alignment_coordinate_state::none>;
    using row_incrementable =
        detail::advanceable_alignment_coordinate<detail::advanceable_alignment_coordinate_state::row>;
    using column_incrementable =
        detail::advanceable_alignment_coordinate<detail::advanceable_alignment_coordinate_state::column>;

    EXPECT_FALSE(std::weakly_incrementable<not_incrementable>);
    EXPECT_TRUE(std::weakly_incrementable<row_incrementable>);
    EXPECT_TRUE(std::weakly_incrementable<column_incrementable>);
}

TEST(advanceable_alignment_coordinate, increment_row)
{
    using row_incrementable =
        detail::advanceable_alignment_coordinate<detail::advanceable_alignment_coordinate_state::row>;

    row_incrementable co{detail::column_index_type{0u}, detail::row_index_type{0u}};
    co = ++co;
    EXPECT_EQ(co.first, 0u);
    EXPECT_EQ(co.second, 1u);
    auto co_tmp = co++;
    EXPECT_EQ(co_tmp.first, 0u);
    EXPECT_EQ(co_tmp.second, 1u);
    EXPECT_EQ(co.first, 0u);
    EXPECT_EQ(co.second, 2u);
    co += 4;
    EXPECT_EQ(co.first, 0u);
    EXPECT_EQ(co.second, 6u);
}

TEST(advanceable_alignment_coordinate, increment_col)
{
    using col_incrementable =
        detail::advanceable_alignment_coordinate<detail::advanceable_alignment_coordinate_state::column>;

    col_incrementable co{detail::column_index_type{0u}, detail::row_index_type{0u}};
    co = ++co;
    EXPECT_EQ(co.first, 1u);
    EXPECT_EQ(co.second, 0u);
    auto co_tmp = co++;
    EXPECT_EQ(co_tmp.first, 1u);
    EXPECT_EQ(co_tmp.second, 0u);
    EXPECT_EQ(co.first, 2u);
    EXPECT_EQ(co.second, 0u);
    co += 4;
    EXPECT_EQ(co.first, 6u);
    EXPECT_EQ(co.second, 0u);
}

TEST(advanceable_alignment_coordinate, decrement_row)
{
    using row_incrementable =
        detail::advanceable_alignment_coordinate<detail::advanceable_alignment_coordinate_state::row>;

    row_incrementable co{detail::column_index_type{0u}, detail::row_index_type{0u}};
    co += 4;
    auto co_tmp = co--;
    EXPECT_EQ(co_tmp.first, 0u);
    EXPECT_EQ(co_tmp.second, 4u);
    EXPECT_EQ(co.first, 0u);
    EXPECT_EQ(co.second, 3u);

    co = --co;
    EXPECT_EQ(co.first, 0u);
    EXPECT_EQ(co.second, 2u);

    co -= 2;
    EXPECT_EQ(co.first, 0u);
    EXPECT_EQ(co.second, 0u);
}

TEST(advanceable_alignment_coordinate, decrement_col)
{
    using col_incrementable =
        detail::advanceable_alignment_coordinate<detail::advanceable_alignment_coordinate_state::column>;

    col_incrementable co{detail::column_index_type{0u}, detail::row_index_type{0u}};
    co += 4;
    auto co_tmp = co--;
    EXPECT_EQ(co_tmp.first, 4u);
    EXPECT_EQ(co_tmp.second, 0u);
    EXPECT_EQ(co.first, 3u);
    EXPECT_EQ(co.second, 0u);

    co = --co;
    EXPECT_EQ(co.first, 2u);
    EXPECT_EQ(co.second, 0u);

    co -= 2;
    EXPECT_EQ(co.first, 0u);
    EXPECT_EQ(co.second, 0u);
}

TEST(advanceable_alignment_coordinate, advance_row)
{
    using row_incrementable =
        detail::advanceable_alignment_coordinate<detail::advanceable_alignment_coordinate_state::row>;

    row_incrementable co{detail::column_index_type{0u}, detail::row_index_type{0u}};

    co = co + 4;
    EXPECT_EQ(co.first, 0u);
    EXPECT_EQ(co.second, 4u);

    co = 4 + co;
    EXPECT_EQ(co.first, 0u);
    EXPECT_EQ(co.second, 8u);
}

TEST(advanceable_alignment_coordinate, advance_col)
{
    using col_incrementable =
        detail::advanceable_alignment_coordinate<detail::advanceable_alignment_coordinate_state::column>;

    col_incrementable co{detail::column_index_type{0u}, detail::row_index_type{0u}};
    co = co + 4;
    EXPECT_EQ(co.first, 4u);
    EXPECT_EQ(co.second, 0u);

    co = 4 + co;
    EXPECT_EQ(co.first, 8u);
    EXPECT_EQ(co.second, 0u);
}

TEST(advanceable_alignment_coordinate, iota_column_index)
{
    using col_incrementable =
        detail::advanceable_alignment_coordinate<detail::advanceable_alignment_coordinate_state::column>;

    col_incrementable co_begin{detail::column_index_type{0u}, detail::row_index_type{0u}};
    col_incrementable co_end{detail::column_index_type{5u}, detail::row_index_type{0u}};
    auto v = std::views::iota(co_begin, co_end);

    EXPECT_TRUE((std::same_as<decltype(v.begin()), decltype(v.end())>));
    EXPECT_EQ((*(--v.end())).first, 4u);

    size_t test = 0u;
    for (auto coordinate : v)
        EXPECT_EQ(coordinate.first, test++);
}

TEST(advanceable_alignment_coordinate, iota_row_index)
{
    using row_incrementable =
        detail::advanceable_alignment_coordinate<detail::advanceable_alignment_coordinate_state::row>;

    row_incrementable co_begin{detail::column_index_type{0u}, detail::row_index_type{0u}};
    row_incrementable co_end{detail::column_index_type{0u}, detail::row_index_type{5u}};
    auto v = std::views::iota(co_begin, co_end);

    EXPECT_TRUE((std::same_as<decltype(v.begin()), decltype(v.end())>));
    EXPECT_EQ((*(--v.end())).second, 4u);

    size_t test = 0u;
    for (auto coordinate : v)
        EXPECT_EQ(coordinate.second, test++);
}

TEST(alignment_coordinate, basic)
{
    EXPECT_TRUE(std::is_default_constructible<alignment_coordinate>::value);
    EXPECT_TRUE(std::is_copy_constructible<alignment_coordinate>::value);
    EXPECT_TRUE(std::is_copy_assignable<alignment_coordinate>::value);
    EXPECT_TRUE(std::is_move_constructible<alignment_coordinate>::value);
    EXPECT_TRUE(std::is_move_assignable<alignment_coordinate>::value);
    EXPECT_TRUE(std::is_destructible<alignment_coordinate>::value);

    using not_incrementable =
        detail::advanceable_alignment_coordinate<detail::advanceable_alignment_coordinate_state::none>;
    using row_incrementable =
        detail::advanceable_alignment_coordinate<detail::advanceable_alignment_coordinate_state::row>;
    using col_incrementable =
        detail::advanceable_alignment_coordinate<detail::advanceable_alignment_coordinate_state::column>;

    not_incrementable co_not{detail::column_index_type{10u}, detail::row_index_type{5u}};
    col_incrementable co_col{detail::column_index_type{10u}, detail::row_index_type{5u}};
    row_incrementable co_row{detail::column_index_type{10u}, detail::row_index_type{5u}};

    alignment_coordinate test1{co_not};
    EXPECT_EQ(test1.first, 10u);
    EXPECT_EQ(test1.second, 5u);

    alignment_coordinate test2{co_col};
    EXPECT_EQ(test2.first, 10u);
    EXPECT_EQ(test2.second, 5u);

    alignment_coordinate test3{co_row};
    EXPECT_EQ(test3.first, 10u);
    EXPECT_EQ(test3.second, 5u);

    alignment_coordinate test4{detail::column_index_type{10u}, detail::row_index_type{5u}};
    EXPECT_EQ(test4.first, 10u);
    EXPECT_EQ(test4.second, 5u);
}

TEST(alignment_coordinate, matrix_coordainte_conversion)
{
    alignment_coordinate co{detail::column_index_type{10u}, detail::row_index_type{5u}};
    detail::matrix_coordinate mc = co;

    EXPECT_EQ(mc.col, 10u);
    EXPECT_EQ(mc.row, 5u);
}

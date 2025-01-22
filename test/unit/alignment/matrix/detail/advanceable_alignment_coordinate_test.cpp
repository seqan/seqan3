// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <type_traits>

#include <seqan3/alignment/matrix/detail/advanceable_alignment_coordinate.hpp>
#include <seqan3/test/expect_same_type.hpp>
#include <seqan3/test/pretty_printing.hpp>

TEST(advanceable_alignment_coordinate, column_index_type)
{
    seqan3::detail::column_index_type ci{1u};
    EXPECT_EQ(ci.get(), 1u);
    EXPECT_SAME_TYPE(std::remove_reference_t<decltype(ci.get())>, size_t);

    seqan3::detail::column_index_type ci2{1};
    EXPECT_EQ(ci2.get(), 1);
    EXPECT_SAME_TYPE(std::remove_reference_t<decltype(ci2.get())>, std::ptrdiff_t);
}

TEST(advanceable_alignment_coordinate, row_index_type)
{
    seqan3::detail::row_index_type ri{1u};
    EXPECT_EQ(ri.get(), 1u);
    EXPECT_SAME_TYPE(std::remove_reference_t<decltype(ri.get())>, size_t);

    seqan3::detail::row_index_type ri2{1};
    EXPECT_EQ(ri2.get(), 1);
    EXPECT_SAME_TYPE(std::remove_reference_t<decltype(ri2.get())>, std::ptrdiff_t);
}

TEST(advanceable_alignment_coordinate, construction)
{
    EXPECT_TRUE(std::is_default_constructible<seqan3::detail::advanceable_alignment_coordinate<>>::value);
    EXPECT_TRUE(std::is_copy_constructible<seqan3::detail::advanceable_alignment_coordinate<>>::value);
    EXPECT_TRUE(std::is_copy_assignable<seqan3::detail::advanceable_alignment_coordinate<>>::value);
    EXPECT_TRUE(std::is_move_constructible<seqan3::detail::advanceable_alignment_coordinate<>>::value);
    EXPECT_TRUE(std::is_move_assignable<seqan3::detail::advanceable_alignment_coordinate<>>::value);
    EXPECT_TRUE(std::is_destructible<seqan3::detail::advanceable_alignment_coordinate<>>::value);
}

TEST(advanceable_alignment_coordinate, construction_with_different_state)
{
    seqan3::detail::advanceable_alignment_coordinate<seqan3::detail::advanceable_alignment_coordinate_state::row> ro{
        seqan3::detail::column_index_type{2u},
        seqan3::detail::row_index_type{3u}};

    seqan3::detail::advanceable_alignment_coordinate<seqan3::detail::advanceable_alignment_coordinate_state::none> no{
        ro};

    EXPECT_EQ(no.first, 2u);
    EXPECT_EQ(no.second, 3u);
}

TEST(advanceable_alignment_coordinate, type_deduction)
{
    using not_incrementable =
        seqan3::detail::advanceable_alignment_coordinate<seqan3::detail::advanceable_alignment_coordinate_state::none>;

    seqan3::detail::advanceable_alignment_coordinate def_co{};
    EXPECT_TRUE((std::is_same_v<decltype(def_co), not_incrementable>));

    seqan3::detail::advanceable_alignment_coordinate co{seqan3::detail::column_index_type{2u},
                                                        seqan3::detail::row_index_type{3u}};
    EXPECT_TRUE((std::is_same_v<decltype(co), not_incrementable>));
}

TEST(advanceable_alignment_coordinate, access)
{
    seqan3::detail::advanceable_alignment_coordinate def_co{};
    EXPECT_EQ(def_co.first, 0u);
    EXPECT_EQ(def_co.second, 0u);

    seqan3::detail::advanceable_alignment_coordinate co{seqan3::detail::column_index_type{2u},
                                                        seqan3::detail::row_index_type{3u}};
    EXPECT_EQ(co.first, 2u);
    EXPECT_EQ(co.second, 3u);
}

TEST(advanceable_alignment_coordinate, weakly_equality_comparable_concept)
{
    using not_incrementable =
        seqan3::detail::advanceable_alignment_coordinate<seqan3::detail::advanceable_alignment_coordinate_state::none>;
    using row_incrementable =
        seqan3::detail::advanceable_alignment_coordinate<seqan3::detail::advanceable_alignment_coordinate_state::row>;
    using column_incrementable = seqan3::detail::advanceable_alignment_coordinate<
        seqan3::detail::advanceable_alignment_coordinate_state::column>;

    EXPECT_TRUE(std::equality_comparable<not_incrementable>);
    EXPECT_TRUE(std::equality_comparable<row_incrementable>);
    EXPECT_TRUE(std::equality_comparable<column_incrementable>);
}

TEST(advanceable_alignment_coordinate, equality)
{
    using test_type =
        seqan3::detail::advanceable_alignment_coordinate<seqan3::detail::advanceable_alignment_coordinate_state::none>;

    test_type t1{seqan3::detail::column_index_type{10u}, seqan3::detail::row_index_type{5u}};
    test_type t2{seqan3::detail::column_index_type{5u}, seqan3::detail::row_index_type{5u}};
    test_type t3{seqan3::detail::column_index_type{10u}, seqan3::detail::row_index_type{10u}};

    EXPECT_TRUE(t1 == t1);
    EXPECT_FALSE(t2 == t1);
    EXPECT_FALSE(t1 == t3);
    EXPECT_FALSE(t2 == t3);
}

TEST(advanceable_alignment_coordinate, inequality)
{
    using test_type =
        seqan3::detail::advanceable_alignment_coordinate<seqan3::detail::advanceable_alignment_coordinate_state::none>;

    test_type t1{seqan3::detail::column_index_type{10u}, seqan3::detail::row_index_type{5u}};
    test_type t2{seqan3::detail::column_index_type{5u}, seqan3::detail::row_index_type{5u}};
    test_type t3{seqan3::detail::column_index_type{10u}, seqan3::detail::row_index_type{10u}};

    EXPECT_FALSE(t1 != t1);
    EXPECT_TRUE(t2 != t1);
    EXPECT_TRUE(t1 != t3);
    EXPECT_TRUE(t2 != t3);
}

TEST(advanceable_alignment_coordinate, incremental_concept)
{
    using not_incrementable =
        seqan3::detail::advanceable_alignment_coordinate<seqan3::detail::advanceable_alignment_coordinate_state::none>;
    using row_incrementable =
        seqan3::detail::advanceable_alignment_coordinate<seqan3::detail::advanceable_alignment_coordinate_state::row>;
    using column_incrementable = seqan3::detail::advanceable_alignment_coordinate<
        seqan3::detail::advanceable_alignment_coordinate_state::column>;

    EXPECT_FALSE(std::weakly_incrementable<not_incrementable>);
    EXPECT_TRUE(std::weakly_incrementable<row_incrementable>);
    EXPECT_TRUE(std::weakly_incrementable<column_incrementable>);
}

TEST(advanceable_alignment_coordinate, increment_row)
{
    using row_incrementable =
        seqan3::detail::advanceable_alignment_coordinate<seqan3::detail::advanceable_alignment_coordinate_state::row>;

    row_incrementable co{seqan3::detail::column_index_type{0u}, seqan3::detail::row_index_type{0u}};
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
    using col_incrementable = seqan3::detail::advanceable_alignment_coordinate<
        seqan3::detail::advanceable_alignment_coordinate_state::column>;

    col_incrementable co{seqan3::detail::column_index_type{0u}, seqan3::detail::row_index_type{0u}};
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
        seqan3::detail::advanceable_alignment_coordinate<seqan3::detail::advanceable_alignment_coordinate_state::row>;

    row_incrementable co{seqan3::detail::column_index_type{0u}, seqan3::detail::row_index_type{0u}};
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
    using col_incrementable = seqan3::detail::advanceable_alignment_coordinate<
        seqan3::detail::advanceable_alignment_coordinate_state::column>;

    col_incrementable co{seqan3::detail::column_index_type{0u}, seqan3::detail::row_index_type{0u}};
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
        seqan3::detail::advanceable_alignment_coordinate<seqan3::detail::advanceable_alignment_coordinate_state::row>;

    row_incrementable co{seqan3::detail::column_index_type{0u}, seqan3::detail::row_index_type{0u}};

    co = co + 4;
    EXPECT_EQ(co.first, 0u);
    EXPECT_EQ(co.second, 4u);

    co = 4 + co;
    EXPECT_EQ(co.first, 0u);
    EXPECT_EQ(co.second, 8u);
}

TEST(advanceable_alignment_coordinate, advance_col)
{
    using col_incrementable = seqan3::detail::advanceable_alignment_coordinate<
        seqan3::detail::advanceable_alignment_coordinate_state::column>;

    col_incrementable co{seqan3::detail::column_index_type{0u}, seqan3::detail::row_index_type{0u}};
    co = co + 4;
    EXPECT_EQ(co.first, 4u);
    EXPECT_EQ(co.second, 0u);

    co = 4 + co;
    EXPECT_EQ(co.first, 8u);
    EXPECT_EQ(co.second, 0u);
}

TEST(advanceable_alignment_coordinate, iota_column_index)
{
    using col_incrementable = seqan3::detail::advanceable_alignment_coordinate<
        seqan3::detail::advanceable_alignment_coordinate_state::column>;

    col_incrementable co_begin{seqan3::detail::column_index_type{0u}, seqan3::detail::row_index_type{0u}};
    col_incrementable co_end{seqan3::detail::column_index_type{5u}, seqan3::detail::row_index_type{0u}};
    auto v = std::views::iota(co_begin, co_end);

    EXPECT_SAME_TYPE(decltype(v.begin()), decltype(v.end()));
    EXPECT_EQ((*(--v.end())).first, 4u);

    size_t test = 0u;
    for (auto coordinate : v)
        EXPECT_EQ(coordinate.first, test++);
}

TEST(advanceable_alignment_coordinate, iota_row_index)
{
    using row_incrementable =
        seqan3::detail::advanceable_alignment_coordinate<seqan3::detail::advanceable_alignment_coordinate_state::row>;

    row_incrementable co_begin{seqan3::detail::column_index_type{0u}, seqan3::detail::row_index_type{0u}};
    row_incrementable co_end{seqan3::detail::column_index_type{0u}, seqan3::detail::row_index_type{5u}};
    auto v = std::views::iota(co_begin, co_end);

    EXPECT_SAME_TYPE(decltype(v.begin()), decltype(v.end()));
    EXPECT_EQ((*(--v.end())).second, 4u);

    size_t test = 0u;
    for (auto coordinate : v)
        EXPECT_EQ(coordinate.second, test++);
}

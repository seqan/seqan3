// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/test/simd_utility.hpp>
#include <seqan3/core/simd/all.hpp>
#include <seqan3/alignment/matrix/alignment_matrix_cell.hpp>

using namespace seqan3;

TEST(affine_matrix_cell_test, scalar_ignore)
{
    seqan3::detail::affine_matrix_cell cell{10, 3, std::ignore};
    auto & [f, s, t] = cell;
    EXPECT_EQ(f, 10);
    EXPECT_EQ(s, 3);
    EXPECT_TRUE(seqan3::detail::decays_to_ignore_v<decltype(t)>);
}

TEST(affine_matrix_cell_test, scalar_trace)
{
    seqan3::detail::affine_matrix_cell cell{10, 3, seqan3::detail::trace_directions::up};
    auto & [f, s, t] = cell;
    EXPECT_EQ(f, 10);
    EXPECT_EQ(s, 3);
    EXPECT_EQ(t, seqan3::detail::trace_directions::up);
}

TEST(affine_matrix_cell_test, simd_ignore)
{
    using native_simd_t = simd_type_t<int32_t>;

    seqan3::detail::affine_matrix_cell cell{simd::fill<native_simd_t>(10), simd::fill<native_simd_t>(3), std::ignore};
    auto & [f, s, t] = cell;
    SIMD_EQ(f, simd::fill<native_simd_t>(10));
    SIMD_EQ(s, simd::fill<native_simd_t>(3));
    EXPECT_TRUE(seqan3::detail::decays_to_ignore_v<decltype(t)>);
}

TEST(affine_matrix_cell_test, simd_trace)
{
    using native_simd_t = simd_type_t<int32_t>;
    using element_type = typename simd_traits<native_simd_t>::scalar_type;

    seqan3::detail::affine_matrix_cell cell{simd::fill<native_simd_t>(10), simd::fill<native_simd_t>(3),
                                            simd::fill<native_simd_t>(
                                                    static_cast<element_type>(detail::trace_directions::up))};
    auto & [f, s, t] = cell;
    SIMD_EQ(f, simd::fill<native_simd_t>(10));
    SIMD_EQ(s, simd::fill<native_simd_t>(3));
    SIMD_EQ(t, simd::fill<native_simd_t>(static_cast<element_type>(detail::trace_directions::up)));
}

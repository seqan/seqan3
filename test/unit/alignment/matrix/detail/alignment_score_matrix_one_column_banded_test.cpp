// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <string>
#include <vector>

#include <seqan3/alignment/matrix/detail/alignment_score_matrix_one_column_banded.hpp>

#include "../../../range/iterator_test_template.hpp"
#include "alignment_matrix_base_test_template.hpp"
#include "simulated_alignment_test_template.hpp"

template <typename t>
struct alignment_score_matrix_one_column_banded_test
{
    using matrix_t = seqan3::detail::alignment_score_matrix_one_column_banded<t>;
    using score_type = t;

    alignment_score_matrix_one_column_banded_test() = default;
    alignment_score_matrix_one_column_banded_test(std::string f, std::string s) :
        matrix{matrix_t{f,
                        s,
                        seqan3::align_cfg::band_fixed_size{seqan3::align_cfg::lower_diagonal{-2},
                                                           seqan3::align_cfg::upper_diagonal{2}},
                        -100}}
    {}

    // Banded matrix. We write only partial columns to the result.
    // The remaining part is therefor 0.
    std::vector<t> gold_matrix{0,  -1, -2, -1, -1, -1, -2, -2, -1, -2, -1, -2, -2,
                               -2, -2, -2, -2, -3, -2, 0,  0,  0,  0,  0,  0};

    matrix_t matrix{};
    size_t last_init_column = 2;
};

INSTANTIATE_TYPED_TEST_SUITE_P(one_column_banded,
                               simulated_alignment_test,
                               alignment_score_matrix_one_column_banded_test<int32_t>, );

using test_type = std::pair<seqan3::detail::alignment_score_matrix_one_column_banded<int32_t>, std::true_type>;

INSTANTIATE_TYPED_TEST_SUITE_P(one_column_banded, alignment_matrix_base_test, test_type, );

//-----------------------------------------------------------------------------
// Test outer iterator
//-----------------------------------------------------------------------------
struct outer_iterator
{};

template <>
struct iterator_fixture<outer_iterator> : alignment_matrix_base_test<test_type>
{
    using base_t = alignment_matrix_base_test<test_type>;

    using iterator_tag = std::input_iterator_tag;
    static constexpr bool const_iterable = false;

    using base_t::test_range;
    std::vector<int> expected_range{0, 0, 0, 0, 0};

    template <typename lhs_t, typename rhs_t>
    static void expect_eq(lhs_t lhs, rhs_t rhs)
    {
        EXPECT_EQ((*lhs.begin()).current, rhs);
    }
};

INSTANTIATE_TYPED_TEST_SUITE_P(banded_score_matrix_outer_iterator, iterator_fixture, outer_iterator, );

//-----------------------------------------------------------------------------
// Test inner iterator
//-----------------------------------------------------------------------------
struct inner_iterator
{};

template <>
struct iterator_fixture<inner_iterator> : alignment_matrix_base_test<test_type>
{
    using base_t = alignment_matrix_base_test<test_type>;

    using iterator_tag = std::forward_iterator_tag;
    static constexpr bool const_iterable = false;

    decltype(*base_t::test_range.begin()) test_range = *base_t::test_range.begin();
    std::vector<int> expected_range{0, 0, 0};

    template <typename lhs_t, typename rhs_t>
    static void expect_eq(lhs_t lhs, rhs_t rhs)
    {
        EXPECT_EQ(lhs.current, rhs);
    }
};

INSTANTIATE_TYPED_TEST_SUITE_P(banded_score_matrix_inner_iterator, iterator_fixture, inner_iterator, );

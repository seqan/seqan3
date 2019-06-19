// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <string>
#include <vector>

#include <seqan3/alignment/matrix/detail/alignment_score_matrix_one_column.hpp>

#include "alignment_score_matrix_test_template.hpp"

using namespace seqan3;

template <typename t>
struct alignment_score_matrix_one_column_test
{
    using matrix_t = detail::alignment_score_matrix_one_column<t>;
    using score_type = t;

    alignment_score_matrix_one_column_test() = default;
    alignment_score_matrix_one_column_test(std::string f, std::string s) : matrix{matrix_t{f, s, -100}}
    {}

    std::vector<t> gold_matrix{  0, -1, -2, -3, -4,
                                -1, -1, -1, -2, -3,
                                -2, -1, -2, -1, -2,
                                -3, -2, -2, -2, -2,
                                -4, -3, -2, -3, -2};

    matrix_t matrix{};
    size_t last_init_column = 4;
};

INSTANTIATE_TYPED_TEST_CASE_P(one_column,
                              alignment_score_matrix_test,
                              alignment_score_matrix_one_column_test<int32_t>);

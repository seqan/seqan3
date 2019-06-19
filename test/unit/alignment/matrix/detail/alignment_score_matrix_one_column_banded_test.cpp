// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <vector>
#include <string>

#include <seqan3/alignment/matrix/detail/alignment_score_matrix_one_column_banded.hpp>

#include "alignment_score_matrix_test_template.hpp"

using namespace seqan3;

template <typename t>
struct alignment_score_matrix_one_column_banded_test
{
    using matrix_t = detail::alignment_score_matrix_one_column_banded<t>;
    using score_type = t;

    alignment_score_matrix_one_column_banded_test() = default;
    alignment_score_matrix_one_column_banded_test(std::string f, std::string s) :
        matrix{matrix_t{f, s, static_band{lower_bound{-2}, upper_bound{2}}, -100}}
    {}

    // Banded matrix. We write only partial columns to the result.
    // The remaining part is therefor 0.
    std::vector<t> gold_matrix{  0, -1, -2,
                                -1, -1, -1, -2,
                                -2, -1, -2, -1, -2,
                                    -2, -2, -2, -2,
                                        -2, -3, -2,
                                 0, 0, 0, 0, 0, 0};

    matrix_t matrix{};
    size_t last_init_column = 2;
};

INSTANTIATE_TYPED_TEST_CASE_P(one_column_banded,
                              alignment_score_matrix_test,
                              alignment_score_matrix_one_column_banded_test<int32_t>);

// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <seqan3/alignment/matrix/edit_distance_trace_matrix_full.hpp>
#include <seqan3/alignment/matrix/row_wise_matrix.hpp>

using score_type = int;
using word_type = uint8_t;

template <bool is_semi_global, bool use_max_errors>
class matrix_type
    : public seqan3::detail::edit_distance_score_matrix_full<word_type, score_type, is_semi_global, use_max_errors>
{
public:
    using base_t = seqan3::detail::edit_distance_score_matrix_full<word_type, score_type,
                                                                   is_semi_global, use_max_errors>;

    matrix_type(size_t const rows_size) : base_t{rows_size}
    {}

    using base_t::max_rows;
    using base_t::add_column;
    using base_t::reserve;
};

static constexpr score_type INF = seqan3::detail::matrix_inf<score_type>;

std::vector<std::vector<int>> as_row_wise_vector(auto matrix)
{
    std::vector<std::vector<int>> result{};
    for (unsigned row = 0; row < matrix.rows(); ++row)
    {
        result.push_back({});
        for (unsigned col = 0; col < matrix.cols(); ++col)
        {
            std::optional<int> entry = matrix.at({seqan3::detail::row_index_type{row},
                                                  seqan3::detail::column_index_type{col}});
            result.back().push_back(entry.value_or(INF));
        }
    }
    return result;
}

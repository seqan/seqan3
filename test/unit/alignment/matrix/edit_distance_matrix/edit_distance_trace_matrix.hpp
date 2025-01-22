// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <seqan3/alignment/matrix/detail/edit_distance_trace_matrix_full.hpp>
#include <seqan3/alignment/matrix/detail/row_wise_matrix.hpp>

#include "../../pairwise/fixture/alignment_fixture.hpp"

using namespace seqan3::test::alignment::fixture; // for N, D, u, l and combinations

using word_type = uint8_t;

template <bool is_semi_global, bool use_max_errors>
class matrix_type : public seqan3::detail::edit_distance_trace_matrix_full<word_type, is_semi_global, use_max_errors>
{
public:
    using base_t = seqan3::detail::edit_distance_trace_matrix_full<word_type, is_semi_global, use_max_errors>;

    matrix_type(size_t const rows_size) : base_t{rows_size}
    {}

    using base_t::add_column;
    using base_t::reserve;
};

std::vector<std::vector<seqan3::detail::trace_directions>> as_row_wise_vector(auto matrix)
{
    std::vector<std::vector<seqan3::detail::trace_directions>> result{};
    for (unsigned row = 0; row < matrix.rows(); ++row)
    {
        result.push_back({});
        for (unsigned col = 0; col < matrix.cols(); ++col)
            result.back().push_back(
                matrix.at({seqan3::detail::row_index_type{row}, seqan3::detail::column_index_type{col}}));
    }
    return result;
}

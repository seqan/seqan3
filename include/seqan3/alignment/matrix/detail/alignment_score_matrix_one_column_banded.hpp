// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::alignment_score_matrix_one_column_banded.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/band/static_band.hpp>
#include <seqan3/alignment/matrix/detail/alignment_matrix_column_major_range_base.hpp>
#include <seqan3/alignment/matrix/detail/alignment_score_matrix_one_column_base.hpp>
#include <seqan3/alignment/matrix/detail/alignment_score_matrix_proxy.hpp>
#include <seqan3/std/iterator>
#include <seqan3/std/ranges>
#include <seqan3/std/span>

namespace seqan3::detail
{

/*!\brief A banded alignment score matrix storing only a single banded column for the computation.
 * \tparam score_t The type of the score; must model either seqan3::arithmetic or seqan3::detail::simd_concept.
 * \ingroup alignment_matrix
 *
 * \details
 *
 * This implementation allocates only a single banded column for the alignment score matrix.
 * The matrix allows access to the underlying values through a range based interface. An iterator over the score matrix
 * iterates in column-major-order. Dereferencing an iterator returns a view over the current matrix column.
 * The value type is seqan3::detail::alignment_score_matrix_proxy, which gives an unified access to the respective
 * matrix cells as needed by the standard alignment algorithm. The matrix is modelled as std::ranges::input_range since
 * the alignment algorithm iterates only once over the complete matrix to calculate the values.
 */
template <typename score_t>
class alignment_score_matrix_one_column_banded :
    protected alignment_score_matrix_one_column_base<score_t>,
    public alignment_matrix_column_major_range_base<alignment_score_matrix_one_column_banded<score_t>>
{
private:
    //!\brief The base class for data storage.
    using matrix_base_t = alignment_score_matrix_one_column_base<score_t>;
    //!\brief The base class for iterating over the matrix.
    using range_base_t = alignment_matrix_column_major_range_base<alignment_score_matrix_one_column_banded<score_t>>;

    //!\brief Befriend the range base class.
    friend range_base_t;

protected:
    using typename matrix_base_t::element_type;
    using typename range_base_t::alignment_column_type;
    //!\copydoc seqan3::detail::alignment_matrix_column_major_range_base::column_data_view_type
    using column_data_view_type = std::span<element_type>;

public:
    using matrix_base_t::num_cols;
    using matrix_base_t::num_rows;

    /*!\name Associated types
     * \{
     */
    //!\copydoc seqan3::detail::alignment_matrix_column_major_range_base::value_type
    using value_type = alignment_score_matrix_proxy<score_t>;
    //!\brief Same as value type.
    using reference = value_type;
    //!\copydoc seqan3::detail::alignment_matrix_column_major_range_base::iterator
    using iterator = typename range_base_t::iterator;
    //!\copydoc seqan3::detail::alignment_matrix_column_major_range_base::sentinel
    using sentinel = typename range_base_t::sentinel;
    using typename matrix_base_t::size_type;
    using typename matrix_base_t::underlying_type;
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    //!\brief Defaulted.
    constexpr alignment_score_matrix_one_column_banded() = default;
    //!\brief Defaulted.
    constexpr alignment_score_matrix_one_column_banded(alignment_score_matrix_one_column_banded const &) = default;
    //!\brief Defaulted.
    constexpr alignment_score_matrix_one_column_banded(alignment_score_matrix_one_column_banded &&) = default;
    //!\brief Defaulted.
    constexpr alignment_score_matrix_one_column_banded &
    operator=(alignment_score_matrix_one_column_banded const &) = default;
    //!\brief Defaulted.
    constexpr alignment_score_matrix_one_column_banded &
    operator=(alignment_score_matrix_one_column_banded &&) = default;
    //!\brief Defaulted.
    ~alignment_score_matrix_one_column_banded() = default;

    /*!\brief Construction from two ranges and a band.
     * \tparam first_sequence_t  The first range type; must model std::ranges::forward_range.
     * \tparam second_sequence_t The second range type; must model std::ranges::forward_range.
     *
     * \param[in] first          The first range.
     * \param[in] second         The second range.
     * \param[in] band           The seqan3::static_band in which to calculate the alignment.
     * \param[in] initial_value  The value to initialise the matrix with. Default initialised if not specified.
     *
     * \details
     *
     * Does only obtain the sizes of the passed ranges in order to allocate the score matrix. Only allocates
     * one column of the size of the band.
     */
    template <std::ranges::forward_range first_sequence_t,
              std::ranges::forward_range second_sequence_t>
    constexpr alignment_score_matrix_one_column_banded(first_sequence_t && first,
                                                       second_sequence_t && second,
                                                       static_band const & band,
                                                       score_t const initial_value = score_t{})
    {
        matrix_base_t::num_cols = static_cast<size_type>(std::ranges::distance(first) + 1);
        matrix_base_t::num_rows = static_cast<size_type>(std::ranges::distance(second) + 1);

        band_col_index = std::min<int32_t>(std::max<int32_t>(band.upper_bound, 0), matrix_base_t::num_cols - 1);
        band_row_index = std::min<int32_t>(std::abs(std::min<int32_t>(band.lower_bound, 0)),
                                           matrix_base_t::num_rows - 1);

        band_size = band_col_index + band_row_index + 1;
        // Reserve one more cell to deal with last cell in the banded column which needs only the diagonal and up cell.
        matrix_base_t::pool.resize(band_size + 1, element_type{initial_value, initial_value});
    }
    //!\}

    //!\brief The column index where the upper bound of the band passes through.
    int32_t band_col_index{};
    //!\brief The row index where the lower bound of the band passes through.
    int32_t band_row_index{};
    //!\brief The size of the band.
    int32_t band_size{};

private:
    //!\copydoc seqan3::detail::alignment_matrix_column_major_range_base::initialise_column
    constexpr alignment_column_type initialise_column(size_type const column_index) noexcept
    {
        int32_t slice_begin = std::max<int32_t>(0, band_col_index - column_index);
        int32_t row_end_index = column_index - band_col_index + band_size;
        int32_t slice_end = band_size - std::max<int32_t>(row_end_index - matrix_base_t::num_rows, 0);

        assert(row_end_index >= 0);
        assert(slice_begin >= 0);
        assert(slice_end > 0);
        assert(slice_begin < slice_end);

        return alignment_column_type{*this,
                                     column_data_view_type{std::addressof(matrix_base_t::pool[slice_begin]),
                                                           std::addressof(matrix_base_t::pool[slice_end])}};
    }

    //!\copydoc seqan3::detail::alignment_matrix_column_major_range_base::make_proxy
    template <std::random_access_iterator iter_t>
    constexpr value_type make_proxy(iter_t host_iter) noexcept
    {
        return {std::get<0>(*host_iter),             // current
                std::get<0>(matrix_base_t::cache),   // last diagonal
                std::get<1>(*(host_iter + 1)),       // last left (read)
                std::get<1>(*(host_iter)),           // next left (write)
                std::get<1>(matrix_base_t::cache)};  // last up
    }

    //!\copydoc seqan3::detail::alignment_matrix_column_major_range_base::on_column_iterator_creation
    template <std::random_access_iterator iter_t>
    constexpr void on_column_iterator_creation(iter_t host_iter) noexcept
    {
        // Cache the last diagonal value.
        std::get<0>(matrix_base_t::cache) = std::get<0>(*host_iter);
    }

    //!\copydoc seqan3::detail::alignment_matrix_column_major_range_base::before_column_iterator_increment
    template <std::random_access_iterator iter_t>
    constexpr void before_column_iterator_increment(iter_t SEQAN3_DOXYGEN_ONLY(host_iter)) noexcept
    {
        // noop.
    }

    //!\copydoc seqan3::detail::alignment_matrix_column_major_range_base::after_column_iterator_increment
    template <std::random_access_iterator iter_t>
    constexpr void after_column_iterator_increment(iter_t host_iter) noexcept
    {
        // Cache the last diagonal value.
        std::get<0>(matrix_base_t::cache) = std::get<0>(*host_iter);
    }
};

} // namespace seqan3::detail

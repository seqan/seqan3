// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::alignment_trace_matrix_full_banded.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/band/static_band.hpp>
#include <seqan3/alignment/matrix/detail/alignment_matrix_column_major_range_base.hpp>
#include <seqan3/alignment/matrix/detail/alignment_trace_matrix_base.hpp>
#include <seqan3/alignment/matrix/detail/alignment_trace_matrix_proxy.hpp>
#include <seqan3/alignment/matrix/detail/trace_iterator_banded.hpp>
#include <seqan3/range/views/zip.hpp>
#include <seqan3/std/iterator>
#include <seqan3/std/ranges>

namespace seqan3::detail
{

/*!\brief An alignment traceback matrix storing the entire banded traceback matrix.
 * \tparam trace_t          The type of the trace directions.
 * \tparam coordinate_only  A boolean flag indicating if only a seqan3::alignment_coordinate should be generated.
 * \ingroup alignment_matrix
 *
 * \details
 *
 * This implementation allocates the full banded traceback matrix using `k*m` memory, where `k` is the band size and
 * `m` the number of columns to store. The matrix allows access to the
 * underlying values through a range based interface. An iterator over the traceback matrix iterates in
 * column-major-order over the traceback matrix. Dereferencing an iterator returns a slice over the current matrix
 * column. The value type is a pair over a seqan3::alignment_coordinate and
 * the seqan3::detail::alignment_trace_matrix_proxy, which gives a unified access to the respective matrix cells
 * as needed by the standard alignment algorithm. The matrix is modelled as std::ranges::input_range since the
 * alignment algorithm iterates only once over the complete matrix to calculate the values.
 *
 * ### Only computing the coordinates
 *
 * Sometimes it is desired to only get access to the alignment coordinates. This can be achieved by setting
 * `coordinate_only = true`. In this case no memory will be allocated and only an internal state is maintained to
 * generate the alignment coordinates. The coordinates are relative to the banded matrix, i.e. they do not represent the
 * corresponding row indices of the actual (unbanded) matrix.
 */
template <typename trace_t, bool coordinate_only = false>
class alignment_trace_matrix_full_banded :
    protected alignment_trace_matrix_base<trace_t>,
    public alignment_matrix_column_major_range_base<alignment_trace_matrix_full_banded<trace_t, coordinate_only>>
{
private:
    static_assert(std::same_as<trace_t, trace_directions> || simd_concept<trace_t>,
                  "Value type must either be a trace_directions object or a simd vector.");

    //!\brief The base class for data storage.
    using matrix_base_t = alignment_trace_matrix_base<trace_t>;
    //!\brief The base class for iterating over the matrix.
    using range_base_t = alignment_matrix_column_major_range_base<alignment_trace_matrix_full_banded<trace_t,
                                                                                                coordinate_only>>;
    //!\brief Befriend the range base class.
    friend range_base_t;

protected:
    using typename matrix_base_t::element_type;
    using typename matrix_base_t::coordinate_type;
    using typename range_base_t::alignment_column_type;
    //!\copydoc alignment_matrix_column_major_range_base::column_data_view_type
    using column_data_view_type = std::conditional_t<coordinate_only,
                                        decltype(std::views::iota(coordinate_type{}, coordinate_type{})),
                                        decltype(views::zip(std::declval<std::span<element_type>>(),
                                                                std::declval<std::span<element_type>>(),
                                                                std::views::iota(coordinate_type{}, coordinate_type{})))>;

public:
    /*!\name Associated types
     * \{
     */
    //!\copydoc seqan3::detail::alignment_matrix_column_major_range_base::value_type
    using value_type = alignment_trace_matrix_proxy<coordinate_type,
                                                    std::conditional_t<coordinate_only,
                                                                       detail::ignore_t const,
                                                                       trace_t>>;
    //!\brief Same as value type.
    using reference = value_type;
    //!\copydoc seqan3::detail::alignment_matrix_column_major_range_base::iterator
    using iterator = typename range_base_t::iterator;
    //!\copydoc seqan3::detail::alignment_matrix_column_major_range_base::sentinel
    using sentinel = typename range_base_t::sentinel;
    using typename matrix_base_t::size_type;
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
     //!\brief Defaulted.
    constexpr alignment_trace_matrix_full_banded() = default;
     //!\brief Defaulted.
    constexpr alignment_trace_matrix_full_banded(alignment_trace_matrix_full_banded const &) = default;
     //!\brief Defaulted.
    constexpr alignment_trace_matrix_full_banded(alignment_trace_matrix_full_banded &&) = default;
     //!\brief Defaulted.
    constexpr alignment_trace_matrix_full_banded & operator=(alignment_trace_matrix_full_banded const &) = default;
     //!\brief Defaulted.
    constexpr alignment_trace_matrix_full_banded & operator=(alignment_trace_matrix_full_banded &&) = default;
     //!\brief Defaulted.
    ~alignment_trace_matrix_full_banded() = default;

    /*!\brief Construction from two ranges and a band.
     * \tparam first_sequence_t  The first range type; must model std::ranges::forward_range.
     * \tparam second_sequence_t The second range type; must model std::ranges::forward_range.
     *
     * \param[in] first         The first range.
     * \param[in] second        The second range.
     * \param[in] band          The seqan3::static_band in which to calculate the alignment.
     * \param[in] initial_value The value to initialise the matrix with. Default initialised if not specified.
     *
     * \details
     *
     * Obtains the sizes of the passed ranges in order to allocate the traceback matrix. Only allocates
     * memory to store the banded matrix.
     * If `coordinate_only` is set to `true`, nothing will be allocated.
     */
    template <std::ranges::forward_range first_sequence_t, std::ranges::forward_range second_sequence_t>
    constexpr alignment_trace_matrix_full_banded(first_sequence_t && first,
                                                 second_sequence_t && second,
                                                 static_band const & band,
                                                 [[maybe_unused]] trace_t const initial_value = trace_t{})
    {
        matrix_base_t::num_cols = static_cast<size_type>(std::ranges::distance(first) + 1);
        matrix_base_t::num_rows = static_cast<size_type>(std::ranges::distance(second) + 1);

        band_col_index = std::min<int32_t>(std::max<int32_t>(band.upper_bound, 0), matrix_base_t::num_cols - 1);
        band_row_index = std::min<int32_t>(std::abs(std::min<int32_t>(band.lower_bound, 0)), matrix_base_t::num_rows - 1);
        band_size = band_col_index + band_row_index + 1;

        // Reserve one more cell to deal with last cell in the banded column which needs only the diagonal and up cell.
        if constexpr (!coordinate_only)
        {
            matrix_base_t::data = typename matrix_base_t::pool_type{number_rows{static_cast<size_type>(band_size)},
                                                                    number_cols{matrix_base_t::num_cols}};
            matrix_base_t::cache_left.resize(band_size + 1, initial_value);
        }
    }
    //!\}

    //!\copydoc seqan3::detail::alignment_trace_matrix_full::trace_path
    auto trace_path(matrix_coordinate const & trace_begin)
    {
        static_assert(!coordinate_only, "Requested trace but storing the trace was disabled!");

        using matrix_iter_t = std::ranges::iterator_t<typename matrix_base_t::pool_type>;
        using trace_iterator_t = trace_iterator_banded<matrix_iter_t>;
        using path_t = std::ranges::subrange<trace_iterator_t, std::ranges::default_sentinel_t>;

        if (trace_begin.row >= static_cast<size_t>(band_size) || trace_begin.col >= matrix_base_t::num_cols)
            throw std::invalid_argument{"The given coordinate exceeds the trace matrix size."};

        return path_t{trace_iterator_t{matrix_base_t::data.begin() + matrix_offset{trace_begin},
                                       column_index_type{band_col_index}},
                      std::ranges::default_sentinel};
    }

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

        coordinate_type row_begin{column_index_type{column_index}, row_index_type{static_cast<size_type>(slice_begin)}};
        coordinate_type row_end{column_index_type{column_index}, row_index_type{static_cast<size_type>(slice_end)}};
        if constexpr (coordinate_only)
        {
            return alignment_column_type{*this,
                                         column_data_view_type{std::views::iota(std::move(row_begin),
                                                                                std::move(row_end))}};
        }
        else
        {
            matrix_coordinate band_begin{row_index_type{static_cast<size_type>(slice_begin)},
                                         column_index_type{column_index}};
            size_type slice_size =  slice_end - slice_begin;
            // We need to jump to the offset.
            auto col = views::zip(
                            std::span<element_type>{std::addressof(matrix_base_t::data[band_begin]), slice_size},
                            std::span<element_type>{std::addressof(matrix_base_t::cache_left[slice_begin]), slice_size},
                            std::views::iota(std::move(row_begin), std::move(row_end)));
            return alignment_column_type{*this, column_data_view_type{std::move(col)}};
        }
    }

    //!\copydoc seqan3::detail::alignment_matrix_column_major_range_base::make_proxy
    template <std::random_access_iterator iter_t>
    constexpr value_type make_proxy(iter_t host_iter) noexcept
    {
        if constexpr (coordinate_only)
        {
            return {*host_iter, std::ignore, std::ignore, std::ignore, std::ignore};
        }
        else
        {
            return {std::get<2>(*host_iter),        // the current coordinate.
                    std::get<0>(*host_iter),        // the current cell.
                    std::get<1>(*(host_iter + 1)),  // the last left cell to read from.
                    std::get<1>(*host_iter),        // the next left cell to write to.
                    matrix_base_t::cache_up         // the last up cell to read/write from/to.
                    };
        }
    }
};

} // namespace seqan3::detail

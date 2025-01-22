// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 * \brief Provides seqan3::detail::edit_distance_trace_matrix_full.
 */

#pragma once

#include <bitset>

#include <seqan3/alignment/matrix/detail/trace_directions.hpp>
#include <seqan3/alignment/pairwise/edit_distance_fwd.hpp>
#include <seqan3/utility/detail/bits_of.hpp>

namespace seqan3::detail
{

/*!\brief The underlying data structure of seqan3::detail::edit_distance_unbanded that represents the
 *        trace matrix.
 * \ingroup alignment_matrix
 * \tparam word_t         \copydoc word_type
 * \tparam is_semi_global \copydoc default_edit_distance_trait_type::is_semi_global
 * \tparam use_max_errors \copydoc default_edit_distance_trait_type::use_max_errors
 */
template <typename word_t, bool is_semi_global, bool use_max_errors>
class edit_distance_trace_matrix_full
{
public:
    //!\brief This friend allows the edit distance algorithm to fill the trace matrix via add_column.
    template <std::ranges::viewable_range database_t,
              std::ranges::viewable_range query_t,
              typename align_config_t,
              typename edit_traits>
    friend class edit_distance_unbanded;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    edit_distance_trace_matrix_full() = default;                                                    //!< Defaulted
    edit_distance_trace_matrix_full(edit_distance_trace_matrix_full const &) = default;             //!< Defaulted
    edit_distance_trace_matrix_full(edit_distance_trace_matrix_full &&) = default;                  //!< Defaulted
    edit_distance_trace_matrix_full & operator=(edit_distance_trace_matrix_full const &) = default; //!< Defaulted
    edit_distance_trace_matrix_full & operator=(edit_distance_trace_matrix_full &&) = default;      //!< Defaulted
    ~edit_distance_trace_matrix_full() = default;                                                   //!< Defaulted

protected:
    //!\brief Allow seqan3::detail::edit_distance_unbanded_trace_matrix_policy to access the private constructor.
    template <typename derived_t, typename edit_traits>
    friend class edit_distance_unbanded_trace_matrix_policy;

    /*!\brief Construct the score_matrix by giving the number of rows within the matrix.
      * \param rows_size \copydoc rows_size
      */
    edit_distance_trace_matrix_full(size_t const rows_size) : rows_size{rows_size}, columns{}
    {}
    //!\}

private:
    struct trace_path_iterator;

public:
    //!\copydoc default_edit_distance_trait_type::word_type
    using word_type = word_t;

    //!\copydoc default_edit_distance_trait_type::word_size
    static constexpr auto word_size = bits_of<word_type>;

    //!\copydoc seqan3::detail::matrix::value_type
    using value_type = detail::trace_directions;

    //!\copydoc seqan3::detail::matrix::reference
    using reference = value_type;

    //!\copydoc seqan3::detail::matrix::size_type
    using size_type = size_t;

    /*!\brief Increase the capacity of the columns to a value that's greater or equal to `new_capacity`.
     * \param new_capacity The new capacity.
     * \details
     *
     * ### Exception
     *
     * Strong exception guarantee.
     */
    void reserve(size_t const new_capacity)
    {
        columns.reserve(new_capacity);
    }

public:
    //!\copydoc seqan3::detail::matrix::at
    reference at(matrix_coordinate const & coordinate) const noexcept
    {
        size_t row = coordinate.row;
        size_t col = coordinate.col;

        assert(row < rows());
        assert(col < cols());

        column_type const & column = columns[col];

        if constexpr (use_max_errors)
            if (!(row < column.max_rows))
                return detail::trace_directions::none;

        if (row == 0u)
        {
            if constexpr (is_semi_global)
                return detail::trace_directions::none;

            if (col == 0u)
                return detail::trace_directions::none;

            return detail::trace_directions::left;
        }

        size_t const idx = (row - 1u) / word_size;
        size_t const offset = (row - 1u) % word_size;

        bool const left = std::bitset<word_size>(column.left[idx])[offset];
        bool const diagonal = std::bitset<word_size>(column.diagonal[idx])[offset];
        bool const up = std::bitset<word_size>(column.up[idx])[offset];

        auto const dir = (left ? detail::trace_directions::left : detail::trace_directions::none)
                       | (diagonal ? detail::trace_directions::diagonal : detail::trace_directions::none)
                       | (up ? detail::trace_directions::up : detail::trace_directions::none);

        return dir;
    }

    //!\copydoc seqan3::detail::matrix::rows
    size_t rows() const noexcept
    {
        return rows_size;
    }

    //!\copydoc seqan3::detail::matrix::cols
    size_t cols() const noexcept
    {
        return columns.size();
    }

    /*!\brief Returns a trace path starting from the given coordinate and ending in the cell with
     *        seqan3::detail::trace_directions::none.
     * \param[in] trace_begin A seqan3::matrix_coordinate pointing to the begin of the trace to follow.
     * \returns A std::ranges::subrange over the corresponding trace path.
     * \throws std::invalid_argument if the specified coordinate is out of range.
     */
    auto trace_path(matrix_coordinate const & trace_begin) const
    {
        if (trace_begin.row >= rows() || trace_begin.col >= cols())
            throw std::invalid_argument{"The given coordinate exceeds the matrix in vertical or horizontal direction."};

        using path_t = std::ranges::subrange<trace_path_iterator, std::default_sentinel_t>;
        return path_t{trace_path_iterator{this, trace_begin}, std::default_sentinel};
    }

protected:
    //!\brief If use_max_errors is true store these additional state information in state_type.
    struct max_errors_state
    {
        //!\brief The number of max_rows within the current column. Computed by
        //!       seqan3::detail::edit_distance_score_matrix_full::max_rows.
        size_t max_rows{};
    };

    //!\brief This information is needed to infer the trace matrix.
    struct trace_matrix_state
    {
        //!\brief Machine words which represent the trace_direction::left.
        std::vector<word_type> left{};
        //!\brief Machine words which represent the trace_direction::diagonal.
        std::vector<word_type> diagonal{};
        //!\brief Machine words which represent the trace_direction::up.
        std::vector<word_type> up{};
    };

    //!\brief The state of one computation step.
    struct column_type : enable_state_t<true, trace_matrix_state>, enable_state_t<use_max_errors, max_errors_state>
    {};

    /*!\brief Adds a column to the trace matrix.
     * \param left     \copydoc trace_matrix_state::left
     * \param diagonal \copydoc trace_matrix_state::diagonal
     * \param up       \copydoc trace_matrix_state::up
     */
    void add_column(std::vector<word_type> left, std::vector<word_type> diagonal, std::vector<word_type> up)
        requires (!use_max_errors)
    {
        column_type column{};
        column.left = std::move(left);
        column.diagonal = std::move(diagonal);
        column.up = std::move(up);

        columns.push_back(std::move(column));
    }

    /*!\brief Adds a column to the trace matrix.
     * \param left     \copydoc trace_matrix_state::left
     * \param diagonal \copydoc trace_matrix_state::diagonal
     * \param up       \copydoc trace_matrix_state::up
     * \param max_rows \copydoc max_errors_state::max_rows
     */
    void add_column(std::vector<word_type> left,
                    std::vector<word_type> diagonal,
                    std::vector<word_type> up,
                    size_t const max_rows)
        requires use_max_errors
    {
        column_type column{};
        column.left = std::move(left);
        column.diagonal = std::move(diagonal);
        column.up = std::move(up);
        column.max_rows = max_rows;

        columns.push_back(std::move(column));
    }

private:
    //!\copydoc seqan3::detail::matrix::rows
    size_t rows_size{};
    //!\brief The columns of the trace matrix.
    std::vector<column_type> columns{};
};

/*!\brief The iterator needed to implement seqan3::detail::edit_distance_trace_matrix_full::trace_path.
 *
 * \details
 *
 * This iterator follows the trace matrix from a starting coordinate until it finds a
 * seqan3::detail::trace_directions::none. This iterator guarantees that it returns exactly only one of
 * these values (normally seqan3::detail::trace_directions can be a combination of these values):
 * * seqan3::detail::trace_directions::left
 * * seqan3::detail::trace_directions::up
 * * seqan3::detail::trace_directions::diagonal
 *
 * This requirement is needed to use the seqan3::detail::aligned_sequence_builder.
 * \extends std::input_iterator
 */
template <typename word_t, bool is_semi_global, bool use_max_errors>
struct edit_distance_trace_matrix_full<word_t, is_semi_global, use_max_errors>::trace_path_iterator
{
    /*!\name Associated types
     * \{
     */
    //!\brief Input iterator tag.
    using iterator_category = std::input_iterator_tag;
    //!\copydoc seqan3::detail::trace_iterator_base::value_type
    using value_type = detail::trace_directions;
    //!\copydoc seqan3::detail::trace_iterator_base::difference_type
    using difference_type = std::ptrdiff_t;
    //!\}

    //!\brief Shortcut for seqan3::detail::trace_directions::diagonal.
    static constexpr value_type D = value_type::diagonal;
    //!\brief Shortcut for seqan3::detail::trace_directions::left.
    static constexpr value_type L = value_type::left;
    //!\brief Shortcut for seqan3::detail::trace_directions::up.
    static constexpr value_type U = value_type::up;
    //!\brief Shortcut for seqan3::detail::trace_directions::none.
    static constexpr value_type N = value_type::none;

    /*!\name Element access
     * \{
     */
    //!\copydoc seqan3::detail::trace_iterator_base::operator*
    constexpr value_type operator*() const
    {
        value_type dir = parent->at(coordinate());

        if (dir == N)
            return N;

        if ((dir & L) == L)
            return L;
        else if ((dir & U) == U)
            return U;
        else
            return D;
    }

    //!\copydoc seqan3::detail::trace_iterator_base::coordinate
    [[nodiscard]] constexpr matrix_coordinate const & coordinate() const
    {
        return coordinate_;
    }
    //!\}

    /*!\name Arithmetic operators
     * \{
     */
    //!\copydoc seqan3::detail::trace_iterator_base::operator++
    constexpr trace_path_iterator & operator++()
    {
        value_type dir = *(*this);

        if ((dir & L) == L)
        {
            coordinate_.col = std::max<size_t>(coordinate_.col, 1) - 1;
        }
        else if ((dir & U) == U)
        {
            coordinate_.row = std::max<size_t>(coordinate_.row, 1) - 1;
        }
        else if ((dir & D) == D)
        {
            coordinate_.row = std::max<size_t>(coordinate_.row, 1) - 1;
            coordinate_.col = std::max<size_t>(coordinate_.col, 1) - 1;
        }

        // seqan3::trace_direction::none in an inner cell of the trace matrix is not possible in
        // non-local alignments, e.g. global, semi-global, max-error, banded. On the other hand,
        // this can happen in a cell of the first row or first colomn.
        assert(dir != N || coordinate_.row == 0 || coordinate_.col == 0);

        return *this;
    }

    //!\copydoc seqan3::detail::trace_iterator_base::operator++
    constexpr void operator++(int)
    {
        ++(*this);
    }
    //!\}

    /*!\name Comparison operators
     * \{
     */
    //!\copydoc seqan3::detail::trace_iterator_base::operator==(derived_t const &, std::default_sentinel_t const &)
    friend bool operator==(trace_path_iterator const & it, std::default_sentinel_t)
    {
        return *it == value_type::none;
    }

    //!\copydoc operator==()
    friend bool operator==(std::default_sentinel_t, trace_path_iterator const & it)
    {
        return it == std::default_sentinel;
    }

    //!\copydoc seqan3::detail::trace_iterator_base::operator!=(derived_t const &, std::default_sentinel_t const &)
    friend bool operator!=(trace_path_iterator const & it, std::default_sentinel_t)
    {
        return !(it == std::default_sentinel);
    }

    //!\copydoc operator!=()
    friend bool operator!=(std::default_sentinel_t, trace_path_iterator const & it)
    {
        return it != std::default_sentinel;
    }
    //!\}

    //!\brief The parent trace matrix.
    edit_distance_trace_matrix_full const * parent{nullptr};
    //!\brief The current coordinate.
    matrix_coordinate coordinate_{};
};

} // namespace seqan3::detail

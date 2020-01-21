// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 * \brief Provides seqan3::detail::edit_distance_trace_matrix_full.
 */

#pragma once

#include <bitset>

#include <seqan3/alignment/matrix/trace_directions.hpp>
#include <seqan3/alignment/pairwise/edit_distance_fwd.hpp>
#include <seqan3/core/bit_manipulation.hpp>

namespace seqan3::detail
{

/*!\brief The underlying data structure of seqan3::detail::edit_distance_unbanded that represents the
 *        trace matrix.
 * \ingroup pairwise_alignment
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
     edit_distance_trace_matrix_full(size_t const rows_size)
         : rows_size{rows_size}, columns{}
     {}
    //!\}

public:
    //!\copydoc default_edit_distance_trait_type::word_type
    using word_type = word_t;

    //!\copydoc default_edit_distance_trait_type::word_size
    static constexpr auto word_size = sizeof_bits<word_type>;

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

        if constexpr(use_max_errors)
            if (!(row < column.max_rows))
                return detail::trace_directions::none;

        if (row == 0u)
        {
            if constexpr(is_semi_global)
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

        auto const dir = (left ? detail::trace_directions::left : detail::trace_directions::none) |
                         (diagonal ? detail::trace_directions::diagonal : detail::trace_directions::none) |
                         (up ? detail::trace_directions::up : detail::trace_directions::none);

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
    struct column_type :
        enable_state_t<true, trace_matrix_state>,
        enable_state_t<use_max_errors, max_errors_state>
    {};

    /*!\brief Adds a column to the trace matrix.
     * \param left     \copydoc trace_matrix_state::left
     * \param diagonal \copydoc trace_matrix_state::diagonal
     * \param up       \copydoc trace_matrix_state::up
     */
    void add_column(std::vector<word_type> left, std::vector<word_type> diagonal, std::vector<word_type> up)
    //!\cond
        requires !use_max_errors
    //!\endcond
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
    void add_column(std::vector<word_type> left, std::vector<word_type> diagonal, std::vector<word_type> up,
                    size_t const max_rows)
    //!\cond
        requires use_max_errors
    //!\endcond
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

} // namespace seqan3::detail

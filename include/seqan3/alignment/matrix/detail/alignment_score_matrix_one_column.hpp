// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::detail::alignment_score_matrix_one_column.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <iterator>
#include <ranges>
#include <span>

#include <seqan3/alignment/matrix/detail/alignment_matrix_column_major_range_base.hpp>
#include <seqan3/alignment/matrix/detail/alignment_score_matrix_one_column_base.hpp>
#include <seqan3/alignment/matrix/detail/alignment_score_matrix_proxy.hpp>
#include <seqan3/utility/concept.hpp>
#include <seqan3/utility/simd/concept.hpp>

namespace seqan3::detail
{

/*!\brief An alignment score matrix storing only a single column for the computation.
 * \tparam score_t The type of the score; must model either seqan3::arithmetic or seqan3::detail::simd_concept.
 * \ingroup alignment_matrix
 *
 * \details
 *
 * This implementation allocates only a single column for the alignment score matrix.
 * The matrix allows access to the underlying values through a range based interface. An iterator over the score matrix
 * iterates in column-major-order. Dereferencing an iterator returns a view over the current matrix column.
 * The value type is seqan3::detail::alignment_score_matrix_proxy, which gives an unified access to the respective
 * matrix cells as needed by the standard alignment algorithm. The matrix is modelled as std::ranges::input_range since
 * the alignment algorithm iterates only once over the complete matrix to calculate the values.
 */
template <typename score_t>
class alignment_score_matrix_one_column :
    protected alignment_score_matrix_one_column_base<score_t>,
    public alignment_matrix_column_major_range_base<alignment_score_matrix_one_column<score_t>>
{
private:
    static_assert(arithmetic<score_t> || simd_concept<score_t>,
                  "The score type must be either an arithmetic type or a simd vector type.");
    //!\brief The base class for data storage.
    using matrix_base_t = alignment_score_matrix_one_column_base<score_t>;
    //!\brief The base class for iterating over the matrix.
    using range_base_t = alignment_matrix_column_major_range_base<alignment_score_matrix_one_column<score_t>>;

    //!\brief Befriend the range base class.
    friend range_base_t;

protected:
    using typename matrix_base_t::element_type;
    using typename range_base_t::alignment_column_type;
    //!\copydoc seqan3::detail::alignment_matrix_column_major_range_base::column_data_view_type
    using column_data_view_type = std::span<element_type>;

    using matrix_base_t::num_cols;

public:
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
    constexpr alignment_score_matrix_one_column() = default;
    //!\brief Defaulted.
    constexpr alignment_score_matrix_one_column(alignment_score_matrix_one_column const &) = default;
    //!\brief Defaulted.
    constexpr alignment_score_matrix_one_column(alignment_score_matrix_one_column &&) = default;
    //!\brief Defaulted.
    constexpr alignment_score_matrix_one_column & operator=(alignment_score_matrix_one_column const &) = default;
    //!\brief Defaulted.
    constexpr alignment_score_matrix_one_column & operator=(alignment_score_matrix_one_column &&) = default;
    //!\brief Defaulted.
    ~alignment_score_matrix_one_column() = default;

    /*!\brief Construction from two ranges.
     * \tparam first_sequence_t  The first range type; must model std::ranges::forward_range.
     * \tparam second_sequence_t The second range type; must model std::ranges::forward_range.
     *
     * \param[in] first         The first range.
     * \param[in] second        The second range.
     * \param[in] initial_value The value to initialise the matrix with. Default initialised if not specified.
     *
     * \details
     *
     * Obtains only the sizes of the passed ranges in order to allocate the score matrix. Only one column
     * is allocated.
     */
    template <std::ranges::forward_range first_sequence_t, std::ranges::forward_range second_sequence_t>
    constexpr alignment_score_matrix_one_column(first_sequence_t && first,
                                                second_sequence_t && second,
                                                score_t const initial_value = score_t{})
    {
        matrix_base_t::num_cols = static_cast<size_type>(std::ranges::distance(first) + 1);
        matrix_base_t::num_rows = static_cast<size_type>(std::ranges::distance(second) + 1);
        matrix_base_t::pool.resize(matrix_base_t::num_rows + 1, element_type{initial_value, initial_value});
    }
    //!\}

private:
    //!\copydoc seqan3::detail::alignment_matrix_column_major_range_base::initialise_column
    constexpr alignment_column_type initialise_column(size_type const SEQAN3_DOXYGEN_ONLY(column_index)) noexcept
    {
        return alignment_column_type{
            *this,
            column_data_view_type{std::addressof(matrix_base_t::pool[0]), matrix_base_t::num_rows}};
    }

    //!\copydoc seqan3::detail::alignment_matrix_column_major_range_base::make_proxy
    template <std::random_access_iterator iter_t>
    constexpr value_type make_proxy(iter_t host_iter) noexcept
    {
        return {std::get<0>(*host_iter),            // current
                std::get<0>(matrix_base_t::cache),  // last diagonal
                std::get<1>(*host_iter),            // last left (read)
                std::get<1>(*host_iter),            // next left (write)
                std::get<2>(matrix_base_t::cache)}; // last up
    }

    //!\copydoc seqan3::detail::alignment_matrix_column_major_range_base::on_column_iterator_creation
    template <std::random_access_iterator iter_t>
    constexpr void on_column_iterator_creation(iter_t host_iter) noexcept
    {
        // Cache the next diagonal value.
        std::get<1>(matrix_base_t::cache) = std::get<0>(*host_iter);
    }

    //!\copydoc seqan3::detail::alignment_matrix_column_major_range_base::before_column_iterator_increment
    template <std::random_access_iterator iter_t>
    constexpr void before_column_iterator_increment(iter_t SEQAN3_DOXYGEN_ONLY(host_iter)) noexcept
    {
        // Update the last diagonal value.
        std::get<0>(matrix_base_t::cache) = std::move(std::get<1>(matrix_base_t::cache));
    }

    //!\copydoc seqan3::detail::alignment_matrix_column_major_range_base::after_column_iterator_increment
    template <std::random_access_iterator iter_t>
    constexpr void after_column_iterator_increment(iter_t host_iter) noexcept
    {
        // Cache the next diagonal value.
        std::get<1>(matrix_base_t::cache) = std::move(std::get<0>(*host_iter));
    }
};

} // namespace seqan3::detail

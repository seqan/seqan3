// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::banded_score_dp_matrix_policy.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <limits>
#include <memory>

#include <range/v3/view/repeat_n.hpp>

#include <seqan3/alignment/matrix/alignment_coordinate.hpp>
#include <seqan3/alignment/matrix/trace_directions.hpp>
#include <seqan3/alignment/pairwise/policy/unbanded_score_dp_matrix_policy.hpp>
#include <seqan3/range/shortcuts.hpp>
#include <seqan3/range/view/slice.hpp>
#include <seqan3/std/ranges>
#include <seqan3/std/span>

namespace seqan3::detail
{

/*!\brief A policy to allocate and manage a banded scoring matrix.
 * \ingroup alingment_policy
 * \tparam derived_t      The type of the derived class.
 * \tparam allocator_type The type of the allocator used to allocate the score matrix.
 */
template <typename derived_t, typename allocator_type>
class banded_score_dp_matrix_policy :
    public unbanded_score_dp_matrix_policy<banded_score_dp_matrix_policy<derived_t, allocator_type>,  allocator_type>
{
private:

    //!\brief The type of the base.
    using base_t = unbanded_score_dp_matrix_policy<banded_score_dp_matrix_policy<derived_t, allocator_type>,  allocator_type>;

    //!\brief Befriend CRTP derived type.
    friend derived_t;

    /*!\name Member types.
     * \{
     */
    //!\brief The type of a matrix cell; inherited from `base_t`.
    using cell_type = typename base_t::cell_type;
    //!\brief The type of the score matrix; inherited from `base_t`.
    using score_matrix_type = typename base_t::score_matrix_type;
    //!\}

    //!\brief A constant value for simulating minus infinity
    static constexpr std::tuple_element_t<0, cell_type> INF =
        std::numeric_limits<std::tuple_element_t<0, cell_type>>::lowest() / 2;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr banded_score_dp_matrix_policy() = default;                                                  //!< Defaulted
    constexpr banded_score_dp_matrix_policy(banded_score_dp_matrix_policy const &) = default;             //!< Defaulted
    constexpr banded_score_dp_matrix_policy(banded_score_dp_matrix_policy &&) = default;                  //!< Defaulted
    constexpr banded_score_dp_matrix_policy & operator=(banded_score_dp_matrix_policy const &) = default; //!< Defaulted
    constexpr banded_score_dp_matrix_policy & operator=(banded_score_dp_matrix_policy &&) = default;      //!< Defaulted
    ~banded_score_dp_matrix_policy() = default;                                                           //!< Defaulted
    //!\}

public:

    /*!\brief Allocates the memory for the dynamic programming matrix given the two sequences.
     * \tparam first_range_t   The type of the first sequence (or packed sequences).
     * \tparam second_range_t  The type of the second sequence (or packed sequences).
     * \tparam band_t          The type of the band.
     * \param[in] first_range  The first sequence (or packed sequences).
     * \param[in] second_range The first sequence (or packed sequences).
     * \param[in] band         The band.
     */
    template <typename first_range_t, typename second_range_t, typename band_t>
    constexpr void allocate_matrix(first_range_t & first_range, second_range_t & second_range, band_t const & band)
    {
        dimension_first_range  = std::ranges::distance(std::ranges::begin(first_range), std::ranges::end(first_range)) + 1;
        dimension_second_range = std::ranges::distance(std::ranges::begin(second_range), std::ranges::end(second_range)) + 1;

        // If upper_bound is negative, set it to 0 and trim the second sequences accordingly.
        band_column_index = std::max(static_cast<uint_fast32_t>(band.upper_bound), static_cast<uint_fast32_t>(0));
        // If lower_bound is positive, set it to 0 and trim the first sequences accordingly.
        band_row_index = std::abs(std::min(static_cast<int_fast32_t>(band.lower_bound),
                                           static_cast<int_fast32_t>(0)));

        // If the band is wider than the sequence length, limit the band width.
        band_column_index = std::min(band_column_index, static_cast<uint_fast32_t>(dimension_second_range - 1));
        band_row_index = std::min(band_row_index, static_cast<uint_fast32_t>(dimension_first_range - 1));

        band_size = band_column_index + band_row_index + 1;

        // Reserve one more cell to deal with last cell in the banded column which needs only the diagonal and up cell.
        // TODO: introduce specific named cell types with initialisation values.
        score_matrix.resize(band_size + 1);

        using std::get;
        get<0>(score_matrix.back()) = INF;
        get<1>(score_matrix.back()) = INF;
        get<2>(score_matrix.back()) = trace_directions::none;
        current_column_index = 0;
        // Position the iterator to the right offset within the band.
        current_matrix_iter = std::ranges::begin(score_matrix) + band_column_index;
    }

    //!\brief Returns the current column of the alignment matrix.
    constexpr auto current_column() noexcept
    {
        auto span = current_band_size();

        assert(span > 0u);  // The span must always be greater than 0.

        advanceable_alignment_coordinate<advanceable_alignment_coordinate_state::row>
            col_begin{column_index_type{current_column_index},
                      row_index_type{static_cast<size_t>(std::ranges::distance(std::ranges::begin(score_matrix),
                                                                               current_matrix_iter))}};
        advanceable_alignment_coordinate<advanceable_alignment_coordinate_state::row>
            col_end{column_index_type{current_column_index}, row_index_type{col_begin.second + span}};

        // Return zip view over current column and current column shifted by one to access the previous horizontal.
        auto zip_score = std::view::zip(std::span{std::addressof(*current_matrix_iter), span},
                                           std::span{std::addressof(*(current_matrix_iter + 1)), span});
        return std::view::zip(std::move(zip_score),
                                 std::view::iota(col_begin, col_end),
                                 ranges::view::repeat_n(std::ignore, span) | std::view::common);
    }

    //!\brief Moves internal matrix pointer to the next column.
    constexpr void go_next_column() noexcept
    {
        // Update the current_column_index.
        base_t::go_next_column();
        // Still in the initialisation phase and need to update the current matrix iterator until begin is reached.
        if (current_matrix_iter != std::ranges::begin(score_matrix))
            --current_matrix_iter;
    }

    //!\brief Returns the current band size depending on the current column position.
    constexpr uint_fast32_t current_band_size() const noexcept
    {
        // Distance from begin of band until end of the entire column (not end of the band).
        int_fast32_t remaining_column_size =
            static_cast<int_fast32_t>(dimension_second_range) -
            std::max(static_cast<int_fast32_t>(0), static_cast<int_fast32_t>(current_column_index - band_column_index));

        // The matrix was trimmed to fit the band exactly, thus the second term cannot be greater or equal
        // than the first term in the equation above.
        assert(remaining_column_size > 0);

        // using const_iter = typename score_matrix_type::const_iterator;
        // The current band size is the min of the remaining column size and the size of the current span of the band.
        return std::min(static_cast<uint_fast32_t>(remaining_column_size),
                        static_cast<uint_fast32_t>(
                                std::ranges::distance(current_matrix_iter, std::ranges::end(score_matrix)) - 1));
    }

    /*!\brief Computes the begin offset of the second_range within the vertical dimension of banded matrix.
     *
     * \details
     *
     * This function can only be called if the current column index is pass the `band_column_index`, i.e. the point
     * where the band does not intersect with the first row of the matrix anymore. Otherwise the result
     * will induce an unsigned underflow and using this return value causes undefined behaviour.
     */
    constexpr uint_fast32_t second_range_begin_offset() const noexcept
    {
        assert(current_column_index > band_column_index);

        return current_column_index - band_column_index - 1;
    }

    //!\brief Checks whether the current band touches the last row.
    constexpr bool band_touches_last_row() const noexcept
    {
        if (current_column_index > band_column_index) // TODO [[likely]]
            return (second_range_begin_offset() + current_band_size() + 1) == dimension_second_range;
        else
            return current_band_size() >= dimension_second_range;
    }

    /*!\brief Trims the sequences to the band parameters.
    * \tparam first_range_t   The type of the first sequence (or packed sequences).
    * \tparam second_range_t  The type of the second sequence (or packed sequences).
    * \tparam band_t          The type of the band.
    * \param[in] first_range  The first sequence (or packed sequences).
    * \param[in] second_range The first sequence (or packed sequences).
    * \param[in] band         The band.
    *
    * \details
    *
    * If the band does not intersect with the origin or the sink of the matrix the
    * sequences are trimmed, such that the band starts in the origin and ends in the sink.
    */
    template <typename first_range_t, typename second_range_t, typename band_t>
    constexpr auto trim_sequences(first_range_t & first_range,
                                  second_range_t & second_range,
                                  band_t const & band) const noexcept
    {
        using band_type = decltype(band.lower_bound);

        band_type dimension_first = std::ranges::distance(std::ranges::begin(first_range), std::ranges::end(first_range));
        band_type dimension_second = std::ranges::distance(std::ranges::begin(second_range), std::ranges::end(second_range));

        auto trim_first_range = [&]() constexpr
        {
            size_t begin_pos = std::max(band.lower_bound - 1, static_cast<band_type>(0));
            size_t end_pos = std::min(band.upper_bound + dimension_second, dimension_first);
            return first_range | view::slice(begin_pos, end_pos);
        };

        auto trim_second_range = [&]() constexpr
        {
            size_t begin_pos = std::abs(std::min(band.upper_bound + 1, static_cast<band_type>(0)));
            size_t end_pos = std::min(dimension_first - band.lower_bound, dimension_second);
            return second_range | view::slice(begin_pos, end_pos);
        };

        return std::tuple{trim_first_range(), trim_second_range()};
    }

    /*!\brief Refines the coordinate for the banded matrix to map the actual sequence position.
     * \param coordinate The coordinate to refine.
     */
    constexpr auto map_banded_coordinate_to_range_position(alignment_coordinate coordinate) const noexcept
    {
        using as_int_t = std::make_signed_t<decltype(coordinate.first)>;
        // Refine the row coordinate to match the original sequence coordinates since the first position of the
        // trace matrix is shifted by the value of the band_column_index, i.e. the upper bound of the band.
        //
        // case 1: ends in column before the band_column_index: subtract the offset from the actual row coordinate.
        // case 2: ends in column after the band_column_index: add the offset to the actual row coordinate.
        coordinate.second += static_cast<as_int_t>(coordinate.first - band_column_index);
        return coordinate;
    }

private:

    using base_t::score_matrix;
    using base_t::dimension_first_range;
    using base_t::dimension_second_range;
    using base_t::current_column_index;

    //!\brief The current matrix iterator.
    typename score_matrix_type::iterator current_matrix_iter;
    //!\brief The column index where the upper bound of the  band starts.
    uint_fast32_t band_column_index{};
    //!\brief The row index where the lower band of the band starts.
    uint_fast32_t band_row_index{};
    //!\brief The full dimension of the band.
    uint_fast32_t band_size{};
};

} // namespace seqan3::detail

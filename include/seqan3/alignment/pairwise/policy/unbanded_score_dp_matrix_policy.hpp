// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::unbanded_score_dp_matrix_policy.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <vector>
#include <tuple>

#include <range/v3/view/iota.hpp>
#include <range/v3/view/repeat_n.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/zip.hpp>

#include <seqan3/alignment/matrix/alignment_coordinate.hpp>
#include <seqan3/core/type_traits/range.hpp>
#include <seqan3/range/shortcuts.hpp>
#include <seqan3/std/span>
#include <seqan3/std/ranges>

namespace seqan3::detail
{

/*!\brief Manages the allocation and provision of an unbanded dynamic programming matrix.
 * \ingroup alignment_policy
 * \tparam derived_t   The derived alignment algorithm.
 * \tparam allocator_t The allocator type used for allocating the dynamic programming matrix.
 */
template <typename derived_t, typename allocator_t>
class unbanded_score_dp_matrix_policy
{
private:

    //!\brief Befriends the derived class to grant it access to the private members.
    friend derived_t;

    /*!\name Member types
     * \{
     */
    //!\brief The underlying cell type of the dynamic programming matrix.
    using cell_type = typename allocator_t::value_type;
    //!\brief The type of the score matrix.
    using score_matrix_type = std::vector<cell_type, allocator_t>;
    //!\}

    /*!\name Constructors, destructor and assignment
     * \brief Defaulted all standard constructor.
     * \{
     */
    constexpr unbanded_score_dp_matrix_policy() = default;                                               //!< Defaulted
    constexpr unbanded_score_dp_matrix_policy(unbanded_score_dp_matrix_policy const &) = default;        //!< Defaulted
    constexpr unbanded_score_dp_matrix_policy(unbanded_score_dp_matrix_policy &&) = default;             //!< Defaulted
    //!\brief Defaulted
    constexpr unbanded_score_dp_matrix_policy & operator=(unbanded_score_dp_matrix_policy const &) = default;
    constexpr unbanded_score_dp_matrix_policy & operator=(unbanded_score_dp_matrix_policy &&) = default; //!< Defaulted
    ~unbanded_score_dp_matrix_policy() = default;                                                        //!< Defaulted
    //!\}

    /*!\brief Allocates the memory for the dynamic programming matrix given the two sequences.
     * \tparam first_range_t   The type of the first sequence (or packed sequences).
     * \tparam second_range_t  The type of the second sequence (or packed sequences).
     * \param[in] first_range  The first sequence (or packed sequences).
     * \param[in] second_range The first sequence (or packed sequences).
     */
    template <typename first_range_t, typename second_range_t>
    constexpr void allocate_matrix(first_range_t & first_range, second_range_t & second_range)
    {
        dimension_first_range = std::ranges::distance(first_range) + 1;
        dimension_second_range = std::ranges::distance(second_range) + 1;

        current_column_index = 0;

        // We use only one column to compute the score.
        score_matrix.resize(dimension_second_range);
    }

    //!\brief Returns the current column of the alignment matrix.
    constexpr auto current_column() noexcept
    {
        advanceable_alignment_coordinate<advanceable_alignment_coordinate_state::row>
            col_begin{column_index_type{current_column_index}, row_index_type{0u}};
        advanceable_alignment_coordinate<advanceable_alignment_coordinate_state::row>
            col_end{column_index_type{current_column_index}, row_index_type{dimension_second_range}};

        return std::view::zip(std::span{score_matrix},
                                 std::view::iota(col_begin, col_end),
                                 ranges::view::repeat_n(std::ignore, dimension_second_range) | std::view::common);
    }

    //!\brief Moves internal matrix pointer to the next column.
    constexpr void go_next_column() noexcept
    {
        ++current_column_index;
    }

    //!\brief The data container.
    score_matrix_type score_matrix{};
    //!\brief Caches the size of the horizontal dimension (number of columns).
    size_t dimension_first_range  = 0;
    //!\brief Caches the size of the vertical dimension (number of rows).
    size_t dimension_second_range = 0;
    //!\brief The index of the active column.
    size_t current_column_index = 0;
};

/*!\brief Returns only the score column of the current matrix column.
 *
 * \details
 *
 * This helper view is used as long as view::get is broken for nested zip-views.
 * See https://github.com/seqan/seqan3/issues/745 for more details.
 */
inline const auto view_get_score_column = ranges::view::transform([](auto && elem)
{
    using std::get;
    return get<0>(std::forward<decltype(elem)>(elem));
});
} // namespace seqan3::detail

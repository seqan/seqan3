// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 * \brief Provides the seqan3::detail::policy_result_builder.
 */

#pragma once

#include <seqan3/search/detail/search_traits.hpp>
#include <seqan3/search/fm_index/concept.hpp>

namespace seqan3::detail
{

//!\brief Provides the function `make_results` if inherited by a search algorithm.
//!\ingroup search
struct policy_result_builder
{
protected:
    /*!\brief Returns all hits (index cursor) without calling locate on each cursor.
     * \tparam index_cursor_t The type of index cursor used in the search algorithm.
     * \tparam configuration_t The search configuration type.
     * \param[in] internal_hits internal_hits A range over internal cursor results.
     * \returns a range over index cursors.
     *
     * \details
     *
     * The result is independent from the search modus (all, single_best, all_best, strata).
     */
    template <typename index_cursor_t, typename configuration_t>
    //!\cond
        requires search_traits<configuration_t>::search_return_index_cursor
    //!\endcond
    auto make_results(std::vector<index_cursor_t> internal_hits, configuration_t const &)
    {
        return internal_hits;
    }

    /*!\brief If `internal_hits` is not empty, calls lazy_locate on the first cursor and returns the first text position.
     * \tparam index_cursor_t The type of index cursor used in the search algorithm.
     * \tparam configuration_t The search configuration type.
     * \param[in] internal_hits internal_hits A range over internal cursor results.
     * \returns a range over a single text position.
     */
    template <typename index_cursor_t, typename configuration_t>
    //!\cond
        requires search_traits<configuration_t>::search_return_text_position &&
                 search_traits<configuration_t>::search_best_hits
    //!\endcond
    auto make_results(std::vector<index_cursor_t> internal_hits, configuration_t const &)
    {
        using index_t = typename index_cursor_t::index_type;
        using position_t = std::conditional_t<index_t::text_layout_mode == text_layout::collection,
                                              std::pair<typename index_t::size_type, typename index_t::size_type>,
                                              typename index_t::size_type>;
        std::vector<position_t> positions;

        if (!internal_hits.empty())
        {
            // only one cursor is reported but it might contain more than one text position
            auto text_pos = internal_hits[0].lazy_locate();
            positions.push_back(text_pos[0]);
        }

        return positions;
    }

    /*!\brief Returns a range over text positions by calling locate on each cursor.
     * \tparam index_cursor_t The type of index cursor used in the search algorithm.
     * \tparam configuration_t The search configuration type.
     * \param[in] internal_hits internal_hits A range over internal cursor results.
     * \returns a range over text positions.
     *
     * \details
     *
     * This function is used for all search modi except single_best (which are all, all_best, and strata).
     *
     * The text positions are sorted and made unique before returning them.
     */
    template <typename index_cursor_t, typename configuration_t>
    //!\cond
        requires search_traits<configuration_t>::search_return_text_position
    //!\endcond
    auto make_results(std::vector<index_cursor_t> internal_hits, configuration_t const &)
    {
        using index_t = typename index_cursor_t::index_type;
        using position_t = std::conditional_t<index_t::text_layout_mode == text_layout::collection,
                                              std::pair<typename index_t::size_type, typename index_t::size_type>,
                                              typename index_t::size_type>;
        std::vector<position_t> positions;

        for (auto const & cursor : internal_hits)
            std::ranges::move(cursor.locate(), std::ranges::back_inserter(positions));

        std::sort(positions.begin(), positions.end());
        positions.erase(std::unique(positions.begin(), positions.end()), positions.end());

        return positions;
    }
};

} // namespace seqan3::detail

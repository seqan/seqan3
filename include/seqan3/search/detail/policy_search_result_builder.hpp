// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 * \brief Provides the seqan3::detail::policy_search_result_builder.
 */

#pragma once

#include <seqan3/core/type_traits/template_inspection.hpp>
#include <seqan3/search/detail/search_traits.hpp>
#include <seqan3/search/fm_index/concept.hpp>

namespace seqan3::detail
{

//!\brief Provides the function `make_results` if inherited by a search algorithm.
//!\ingroup search
template <typename search_configuration_t>
//!\cond
    requires is_type_specialisation_of_v<search_configuration_t, configuration>
//!\endcond
struct policy_search_result_builder
{
protected:
    //!\brief The traits type over the search configuration.
    using search_traits_type = detail::search_traits<search_configuration_t>;
    //!\brief The configured search result type.
    using callback_type = typename search_traits_type::search_result_type;
    //!\brief The actual search result type.
    using search_result_type = std::tuple_element_t<1, callback_type>;

    static_assert(!std::same_as<callback_type, empty_type>, "The search result type was not configured properly.");

    /*!\brief Returns all hits (index cursor) without calling locate on each cursor.
     * \tparam index_cursor_t The type of index cursor used in the search algorithm.
     * \param[in] internal_hits internal_hits A range over internal cursor results.
     * \returns a range over index cursors.
     *
     * \details
     *
     * The result is independent from the search modus (all, single_best, all_best, strata).
     */
    template <typename index_cursor_t>
    auto make_results(std::vector<index_cursor_t> internal_hits)
    {
        return internal_hits;
    }

    /*!\brief If `internal_hits` is not empty, calls lazy_locate on the first cursor and returns the first text position.
     * \tparam index_cursor_t The type of index cursor used in the search algorithm.
     * \param[in] internal_hits internal_hits A range over internal cursor results.
     * \returns a range over a single text position.
     */
    template <typename index_cursor_t>
    //!\cond
        requires search_traits_type::search_return_text_position &&
                 search_traits_type::search_single_best_hit
    //!\endcond
    auto make_results(std::vector<index_cursor_t> internal_hits)
    {
        std::vector<search_result_type> positions{};

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
     * \param[in] internal_hits internal_hits A range over internal cursor results.
     * \returns a range over text positions.
     *
     * \details
     *
     * This function is used for all search modi except single_best (which are all, all_best, and strata).
     *
     * The text positions are sorted and made unique before returning them.
     */
    template <typename index_cursor_t>
    //!\cond
        requires search_traits_type::search_return_text_position &&
                 (!search_traits_type::search_single_best_hit)
    //!\endcond
    auto make_results(std::vector<index_cursor_t> internal_hits)
    {
        std::vector<search_result_type> positions{};

        for (auto const & cursor : internal_hits)
            std::ranges::move(cursor.locate(), std::ranges::back_inserter(positions));

        std::sort(positions.begin(), positions.end());
        positions.erase(std::unique(positions.begin(), positions.end()), positions.end());

        return positions;
    }
};

} // namespace seqan3::detail

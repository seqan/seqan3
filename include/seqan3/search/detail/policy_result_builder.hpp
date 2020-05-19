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
#include <seqan3/search/search_result.hpp>

namespace seqan3::detail
{

//!\brief Provides the function `make_results` if inherited by a search algorithm.
//!\ingroup search
struct policy_result_builder
{
protected:
    /*!\brief Returns all hits (index cursors) without calling locate on each cursor.
     * \tparam index_cursor_t The type of index cursor used in the search algorithm.
     * \tparam configuration_t The search configuration type.
     * \param[in] internal_hits internal_hits A range over internal cursor results.
     * \returns a range over seqan3::search_result's.
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
        using search_result_t = search_result<size_t, index_cursor_t, size_t, size_t>;
        std::vector<search_result_t> results;
        results.reserve(internal_hits.size());

        for (auto const & cursor : internal_hits)
            results.push_back(search_result_t{0, cursor});

        return results;
    }

    /*!\brief If `internal_hits` is not empty, calls lazy_locate on the first cursor and returns a
    *         seqan3::search_result with the first text position.
     * \tparam index_cursor_t The type of index cursor used in the search algorithm.
     * \tparam configuration_t The search configuration type.
     * \param[in] internal_hits internal_hits A range over internal cursor results.
     * \returns a range over a seqan3::search_result.
     */
    template <typename index_cursor_t, typename configuration_t>
    //!\cond
        requires search_traits<configuration_t>::search_return_text_position &&
                 search_traits<configuration_t>::search_single_best_hit
    //!\endcond
    auto make_results(std::vector<index_cursor_t> internal_hits, configuration_t const &)
    {
        using index_size_t = typename index_cursor_t::index_type::size_type;
        using search_result_t = search_result<size_t, empty_type, size_t, index_size_t>;
        std::vector<search_result_t> results;

        if (!internal_hits.empty())
        {
            // only one cursor is reported but it might contain more than one text position
            auto [ref_id, ref_pos] = internal_hits[0].lazy_locate()[0];
            results.push_back(search_result_t{0, ref_id, ref_pos});
        }
        return results;
    }

    /*!\brief Returns a range over seqan3::search_result's by calling locate on each cursor.
     * \tparam index_cursor_t The type of index cursor used in the search algorithm.
     * \tparam configuration_t The search configuration type.
     * \param[in] internal_hits internal_hits A range over internal cursor results.
     * \returns a range over seqan3::search_result's.
     *
     * \details
     *
     * This function is used for all search modi except single_best (which are all, all_best, and strata).
     *
     * The text positions are sorted and made unique by position before returning them.
     */
    template <typename index_cursor_t, typename configuration_t>
    //!\cond
        requires search_traits<configuration_t>::search_return_text_position &&
                 (!search_traits<configuration_t>::search_single_best_hit)
    //!\endcond
    auto make_results(std::vector<index_cursor_t> internal_hits, configuration_t const &)
    {
        using index_size_t = typename index_cursor_t::index_type::size_type;
        using search_result_t = search_result<size_t, empty_type, size_t, index_size_t>;

        std::vector<search_result_t> results;
        results.reserve(internal_hits.size()); // we at least as many text positions as cursors, possibly more

        for (auto const & cursor : internal_hits)
        {
            for (auto [ref_id, ref_pos] : cursor.locate())
            {
                results.push_back(search_result_t{0, ref_id, ref_pos});
                auto res = results.back();
            }
        }

        // sort by reference id and reference position
        auto compare = [] (auto const & r1, auto const & r2)
        {
            if (r1.reference_id() == r2.reference_id())
                return r1.reference_begin_pos() < r2.reference_begin_pos();
            return r1.reference_id() < r2.reference_id();
        };

        std::sort(results.begin(), results.end(), compare);
        results.erase(std::unique(results.begin(), results.end()), results.end());

        return results;
    }
};

} // namespace seqan3::detail

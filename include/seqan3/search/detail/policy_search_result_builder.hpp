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
#include <seqan3/search/search_result.hpp>

namespace seqan3::detail
{

//!\brief Provides the function `make_results` if inherited by a search algorithm.
//!\ingroup search
template <typename search_configuration_t>
#if !SEQAN3_WORKAROUND_GCC_93467
//!\cond
    requires is_type_specialisation_of_v<search_configuration_t, configuration>
//!\endcond
#endif // !SEQAN3_WORKAROUND_GCC_93467
struct policy_search_result_builder
{
protected:
    //!\brief The traits type over the search configuration.
    using search_traits_type = detail::search_traits<search_configuration_t>;
    //!\brief The configured search result type.
    using search_result_type = typename search_traits_type::search_result_type;

    static_assert(!std::same_as<search_result_type, empty_type>, "The search result type was not configured properly.");

    /*!\brief Returns all hits (index cursors) without calling locate on each cursor.
     *
     * \tparam index_cursor_t The type of index cursor used in the search algorithm.
     * \tparam query_index_t The index type of the query.
     * \tparam callback_t The callback which is called for every hit.
     *
     * \param[in] internal_hits internal_hits A range over internal cursor results.
     * \param[in] idx The index associated with the current query.
     * \param[in] callback The callback to invoke for every hit.
     *
     * \details
     *
     * The result is independent from the search modus (all, single_best, all_best, strata).
     */
    template <typename index_cursor_t, typename query_index_t, typename callback_t>
    void make_results(std::vector<index_cursor_t> internal_hits, query_index_t idx, callback_t && callback)
    {
        static_assert(std::destructible<decltype(search_result_type{std::declval<query_index_t &&>(),
                                                                    std::declval<index_cursor_t &>()})>,
                      "The configured search result type is not compatible with the requested search output. "
                      "Note trying to instantiate search result with an index cursor.");

        for (size_t i = 0; i < internal_hits.size(); ++i)
            callback(search_result_type{std::move(idx), internal_hits[i]});
    }

    /*!\brief If `internal_hits` is not empty, calls lazy_locate on the first cursor and returns a
    *         seqan3::search_result with the first text position.
    *
    * \tparam index_cursor_t The type of index cursor used in the search algorithm.
    * \tparam query_index_t The index type of the query.
    * \tparam callback_t The callback which is called for every hit.
    *
    * \param[in] internal_hits internal_hits A range over internal cursor results.
    * \param[in] idx The index associated with the current query.
    * \param[in] callback The callback to invoke for every hit.
    */
    template <typename index_cursor_t, typename query_index_t, typename callback_t>
    //!\cond
        requires search_traits_type::search_return_text_position &&
                 search_traits_type::search_single_best_hit
    //!\endcond
    void make_results(std::vector<index_cursor_t> internal_hits, query_index_t idx, callback_t && callback)
    {
        using ref_location_pair_t = decltype(std::declval<index_cursor_t &>().lazy_locate().front());
        static_assert(std::destructible<decltype(search_result_type{std::declval<query_index_t &&>(),
                                                                    std::declval<ref_location_pair_t &>().first,
                                                                    std::declval<ref_location_pair_t &>().second})>,
                      "The configured search result type is not compatible with the requested search output. "
                      "Note trying to instantiate search result with text positions and single best hit strategy.");

        if (!internal_hits.empty())
        {
            // only one cursor is reported but it might contain more than one text position
            auto && [ref_id, ref_pos] = internal_hits[0].lazy_locate()[0];
            callback(search_result_type{std::move(idx), ref_id, ref_pos});
        }
    }

    /*!\brief Returns a range over seqan3::search_result by calling locate on each cursor.
     *
     * \tparam index_cursor_t The type of index cursor used in the search algorithm.
     * \tparam query_index_t The index type of the query.
     * \tparam callback_t The callback which is called for every hit.
     *
     * \param[in] internal_hits internal_hits A range over internal cursor results.
     * \param[in] idx The index associated with the current query.
     * \param[in] callback The callback to invoke for every hit.
     *
     * \details
     *
     * This function is used for all search modi except single_best (which are all, all_best, and strata).
     *
     * The text positions are sorted and made unique by position before invoking the callback on them.
     */
    template <typename index_cursor_t, typename query_index_t, typename callback_t>
    //!\cond
        requires search_traits_type::search_return_text_position &&
                 (!search_traits_type::search_single_best_hit)
    //!\endcond
    void make_results(std::vector<index_cursor_t> internal_hits, query_index_t idx, callback_t && callback)
    {
        using ref_location_pair_t = decltype(std::declval<index_cursor_t>().lazy_locate().front());
        static_assert(std::destructible<decltype(search_result_type{std::declval<query_index_t &&>(),
                                                                    std::declval<ref_location_pair_t &>().first,
                                                                    std::declval<ref_location_pair_t &>().second})>,
                      "The configured search result type is not compatible with the requested search output. "
                      "Note trying to instantiate search result with text positions with either all, all_best or "
                      "strata hit strategy.");

        std::vector<search_result_type> results{};
        results.reserve(internal_hits.size()); // expect at least as many text positions as cursors, possibly more

        for (auto const & cursor : internal_hits)
            for (auto && [ref_id, ref_pos] : cursor.locate())
                results.push_back(search_result_type{std::move(idx), ref_id, ref_pos});

        // sort by reference id or by reference position if both have the same reference id.
        auto compare = [] (auto const & r1, auto const & r2)
        {
            return (r1.reference_id() == r2.reference_id()) ? (r1.reference_begin_pos() < r2.reference_begin_pos())
                                                            : (r1.reference_id() < r2.reference_id());
        };

        std::sort(results.begin(), results.end(), compare);
        results.erase(std::unique(results.begin(), results.end()), results.end());

        for (auto && search_result : results)
            callback(std::move(search_result));
    }
};

} // namespace seqan3::detail

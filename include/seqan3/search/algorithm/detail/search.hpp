// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
 * \brief Provides the public interface for search algorithms.
 */

#pragma once

#include <seqan3/core/type_traits/pre.hpp>
#include <seqan3/search/algorithm/detail/search_scheme_algorithm.hpp>
#include <seqan3/search/algorithm/detail/search_traits.hpp>
#include <seqan3/search/algorithm/detail/search_trivial.hpp>
#include <seqan3/search/configuration/all.hpp>
#include <seqan3/search/fm_index/concept.hpp>

namespace seqan3::detail
{

/*!\addtogroup submodule_search_algorithm
 * \{
 */

/*!\brief Search a single query in an index.
 * \tparam index_t   Must model seqan3::fm_index_specialisation.
 * \tparam queries_t Must model std::ranges::random_access_range over the index's alphabet.
 * \param[in] index  String index to be searched.
 * \param[in] query  A single query.
 * \param[in] cfg    A configuration object specifying the search parameters.
 * \returns `True` if and only if `abort_on_hit` is `true` and a hit has been found.
 *
 * ### Complexity
 *
 * \f$O(|query|^e)\f$ where \f$e\f$ is the maximum number of errors.
 *
 * ### Exceptions
 *
 * Strong exception guarantee if iterating the query does not change its state and if invoking a possible delegate
 * specified in `cfg` also has a strong exception guarantee; basic exception guarantee otherwise.
 */
template <typename index_t, typename query_t, typename configuration_t, typename hits_t>
// TODO add another separate argument for cursors
inline auto search_single(index_t const & index, query_t & query, configuration_t const & cfg, hits_t & hits)
{
    using search_traits_t = search_traits<configuration_t>;

    // retrieve error numbers / rates
    detail::search_param max_error{0, 0, 0, 0};
    if constexpr (search_traits_t::search_with_max_error)
    {
        auto & [total, subs, ins, del] = max_error;
        std::tie(total, subs, ins, del) = std::apply([](auto ...args){ return std::tuple{args...}; },
                                                     get<search_cfg::max_error>(cfg).value);
    }
    else if constexpr (search_traits_t::search_with_max_error_rate)
    {
        // NOTE: Casting doubles rounds towards zero (i.e. floor for positive numbers). Thus, given a rate of
        // 10% and a read length of 101 the max number of errors is correctly casted from 10.1 errors to 10
        auto & [total, subs, ins, del] = max_error;
        std::tie(total, subs, ins, del) = std::apply([& query] (auto && ... args)
            {
                return std::tuple{(args * std::ranges::size(query))...};
            }, get<search_cfg::max_error_rate>(cfg).value);
    }

    // TODO: if total not set: max_error.total = max_error.deletion + max_error.substitution + max_error.insertion;
    // TODO: throw exception when any error number or rate is higher than the total error number/rate
    // throw std::invalid_argument("The total number of errors is set to zero while there is a positive number"
    //                             " of errors for a specific error type.");

    // construct internal delegate for collecting hits for later filtering (if necessary)

    //TODO this needs to be externalised as well, but logic in here needs to change a little (never clear in here!)
    std::vector<typename index_t::cursor_type> internal_hits;
    auto internal_delegate = [&internal_hits, &max_error] (auto const & it)
    {
        internal_hits.push_back(it);
    };

    // choose mode
    if constexpr (search_traits_t::search_best_hits)
    {
        detail::search_param max_error2{max_error};
        max_error2.total = 0;
        while (internal_hits.empty() && max_error2.total <= max_error.total)
        {
            detail::search_algo<true>(index, query, max_error2, internal_delegate);
            max_error2.total++;
        }
    }
    else if constexpr (search_traits_t::search_all_best_hits)
    {
        detail::search_param max_error2{max_error};
        max_error2.total = 0;
        while (internal_hits.empty() && max_error2.total <= max_error.total)
        {
            detail::search_algo<false>(index, query, max_error2, internal_delegate);
            max_error2.total++;
        }
    }
    else if constexpr (search_traits_t::search_strata_hits)
    {
        detail::search_param max_error2{max_error};
        max_error2.total = 0;
        while (internal_hits.empty() && max_error2.total <= max_error.total)
        {
            detail::search_algo<true>(index, query, max_error2, internal_delegate);
            max_error2.total++;
        }
        if (!internal_hits.empty())
        {
            internal_hits.clear(); // TODO: don't clear when using Optimum Search Schemes with lower error bounds
            uint8_t const s = get<search_cfg::mode>(cfg).value;
            max_error2.total += s - 1;
            detail::search_algo<false>(index, query, max_error2, internal_delegate);
        }
    }
    else // detail::search_mode_all
    {
        detail::search_algo<false>(index, query, max_error, internal_delegate);
    }

    // TODO: filter hits and only do it when necessary (depending on error types)

    // output cursors or text_positions
    if constexpr (search_traits_t::search_return_index_cursor)
    {
        hits = std::move(internal_hits);
    }
    else
    {
        if constexpr (search_traits_t::search_best_hits)
        {

            // only one cursor is reported but it might contain more than one text position
            if (!internal_hits.empty())
            {
                auto text_pos = internal_hits[0].lazy_locate();
                hits.push_back(text_pos[0]);
            }
        }
        else
        {
            hits.reserve(internal_hits.size());
            for (auto const & cur : internal_hits)
                cur.locate(hits);

            //TODO when are hits not unique?
            if (max_error.insertion || max_error.deletion)
            {
                std::ranges::sort(hits);
                hits.erase(std::unique(hits.begin(), hits.end()), hits.end());
            }
        }
    }
}

/*!\brief Search a query or a range of queries in an index.
 * \tparam index_t    Must model seqan3::fm_index_specialisation.
 * \tparam queries_t  Must model std::ranges::random_access_range over the index's alphabet.
 *                    a range of queries must additionally model std::ranges::forward_range.
 * \param[in] index   String index to be searched.
 * \param[in] queries A single query or a range of queries.
 * \param[in] cfg     A configuration object specifying the search parameters.
 * \returns `True` if and only if `abort_on_hit` is `true` and a hit has been found.
 *
 * ### Complexity
 *
 * Each query takes \f$O(|query|^e)\f$ where \f$e\f$ is the maximum number of errors.
 *
 * ### Exceptions
 *
 * Strong exception guarantee if iterating the query does not change its state and if invoking a possible delegate
 * specified in `cfg` also has a strong exception guarantee; basic exception guarantee otherwise.
 */
template <typename index_t, typename queries_t, typename configuration_t>
inline auto search_all(index_t const & index, queries_t & queries, configuration_t const & cfg)
{
    using cfg_t = remove_cvref_t<configuration_t>;
    // return type: for each query: a vector of text_positions (or cursors)
    // delegate params: text_position (or cursor). we will withhold all hits of one query anyway to filter
    //                  duplicates. more efficient to call delegate once with one vector instead of calling
    //                  delegate for each hit separately at once.
    using text_pos_t = std::conditional_t<index_t::text_layout_mode == text_layout::collection,
                                          std::pair<typename index_t::size_type, typename index_t::size_type>,
                                          typename index_t::size_type>;
    using hit_t = std::conditional_t<cfg_t::template exists<search_cfg::output<detail::search_output_index_cursor>>(),
                                     typename index_t::cursor_type,
                                     text_pos_t>;

    //TODO have separate object here for cursors and for hits
    // make these externalisable via a config object; config object should hold thread_local variables to be thread-safe
    if constexpr (std::ranges::forward_range<queries_t> &&
                  std::ranges::random_access_range<std::ranges::range_reference_t<queries_t>>)
    {
        // TODO: if constexpr (contains<search_cfg::id::on_hit>(cfg))
        std::vector<std::vector<hit_t>> hits;
        hits.resize(std::ranges::distance(queries));
        size_t count = 0;
        for (auto && query : queries)
        {
            search_single(index, query, cfg, hits[count]);
            ++count;
        }
        return hits;
    }
    else // std::ranges::random_access_range<queries_t>
    {
        std::vector<hit_t> hits;
        // TODO: if constexpr (contains<search_cfg::id::on_hit>(cfg))
        search_single(index, queries, cfg, hits);
        return hits;
    }
}

//!\}

} // namespace seqan3::detail

// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
 * \brief Provides the public interface for search algorithms.
 */

#pragma once

#include <seqan3/core/type_traits/pre.hpp>
#include <seqan3/search/algorithm/configuration/all.hpp>
#include <seqan3/search/algorithm/detail/search_scheme_algorithm.hpp>
#include <seqan3/search/algorithm/detail/search_trivial.hpp>
#include <seqan3/search/fm_index/concept.hpp>

namespace seqan3::detail
{

/*!\addtogroup submodule_search_algorithm
 * \{
 */

/*!\brief Search a single query in an index.
 * \tparam index_t   Must model seqan3::FmIndex.
 * \tparam queries_t Must model std::ranges::RandomAccessRange over the index's alphabet.
 * \param[in] index  String index to be searched.
 * \param[in] text   String to be searched. (Optional: Can be empty vector or collection.)
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
template <typename index_t, typename text_t, typename query_t, typename configuration_t>
inline auto search_single(index_t const & index, text_t const & text, query_t & query, configuration_t const & cfg)
{
    using cfg_t = remove_cvref_t<configuration_t>;

    // retrieve error numbers / rates
    detail::search_error_param max_error{0, 0, 0, 0};
    if constexpr (cfg.template exists<search_cfg::max_error>())
    {

        auto & [total, subs, ins, del] = max_error;
        std::tie(total, subs, ins, del) = std::apply([](auto ...args){ return std::tuple{args...}; },
                                                     get<search_cfg::max_error>(cfg).value);
    }
    else if constexpr (cfg.template exists<search_cfg::max_error_rate>())
    {
        // NOTE: Casting doubles rounds towards zero (i.e. floor for positive numbers). Thus, given a rate of
        // 10% and a read length of 101 the max number of errors is correctly casted from 10.1 errors to 10
        auto & [total, subs, ins, del] = max_error;
        std::tie(total, subs, ins, del) = std::apply([& query] (auto && ... args)
            {
                return std::tuple{(args * std::ranges::size(query))...};
            }, get<search_cfg::max_error_rate>(cfg).value);
    }

    detail::search_param search_params{max_error};

    if constexpr (cfg.template exists<search_cfg::itv_threshold>())
    {
        search_params.itv_threshold = get<search_cfg::itv_threshold>(cfg).value.first;
        search_params.min_step = get<search_cfg::itv_threshold>(cfg).value.second;
    }
    else
    {
        //TODO add some formula with text length alphabet size to calculate default itv_threshold
        search_params.itv_threshold = 10;

        //TODO test for aminoacids
        if constexpr (dimension_v<text_t> != 1){
            uint64_t cumlength = 0;
            for(size_t i = 0; i < text.size(); ++i)
                cumlength += text[i].size();
            search_params.min_step = std::round(log(cumlength) / log(4)) + 4; //+ 3 or 4
        }
        else
        {
            search_params.min_step = std::round(log(text.size()) / log(4)) + 4;
        }
    }

    // TODO: if total not set: max_error.total = max_error.deletion + max_error.substitution + max_error.insertion;
    // TODO: throw exception when any error number or rate is higher than the total error number/rate
    // throw std::invalid_argument("The total number of errors is set to zero while there is a positive number"
    //                             " of errors for a specific error type.");

    // construct internal delegate for collecting hits for later filtering (if necessary)
    std::vector<typename index_t::cursor_type> internal_hits;
    auto internal_delegate = [&internal_hits, &search_params] (auto const & it)
    {
        internal_hits.push_back(it);
    };

     using hit_t = std::conditional_t<index_t::is_collection_,
                                         std::pair<typename index_t::size_type, typename index_t::size_type>,
                                         typename index_t::size_type>;

    std::vector<hit_t> internal_itv_hits;
    auto internal_itv_delegate = [&internal_itv_hits, &search_params] (auto const & pos)
    {
        internal_itv_hits.push_back(pos);
    };

    // choose mode
    if constexpr (cfg_t::template exists<search_cfg::mode<detail::search_mode_best>>())
    {
        detail::search_param search_params2{search_params};
        search_params2.total = 0;
        while (internal_hits.empty() && search_params2.total <= search_params.total)
        {
            detail::search_algo<true>(index, text, query, search_params2, internal_delegate, internal_itv_delegate);
            search_params2.total++;
        }
    }
    else if constexpr (cfg_t::template exists<search_cfg::mode<detail::search_mode_all_best>>())
    {
        detail::search_param search_params2{search_params};
        search_params2.total = 0;
        while (internal_hits.empty() && search_params2.total <= search_params.total)
        {
            detail::search_algo<false>(index, text, query, search_params2, internal_delegate, internal_itv_delegate);
            search_params2.total++;
        }
    }
    else if constexpr (cfg_t::template exists<search_cfg::mode<search_cfg::strata>>())
    {
        detail::search_param search_params2{search_params};
        search_params2.total = 0;
        while (internal_hits.empty() && search_params2.total <= search_params.total)
        {
            detail::search_algo<true>(index, text, query, search_params2, internal_delegate, internal_itv_delegate);
            search_params2.total++;
        }
        if (!internal_hits.empty())
        {
            internal_hits.clear(); // TODO: don't clear when using Optimum Search Schemes with lower error bounds
            uint8_t const s = get<search_cfg::mode>(cfg).value;
            search_params2.total += s - 1;
            detail::search_algo<false>(index, text, query, search_params2, internal_delegate, internal_itv_delegate);
        }
    }
    else // detail::search_mode_all
    {
        detail::search_algo<false>(index, text, query, search_params, internal_delegate, internal_itv_delegate);
    }

    // TODO: filter hits and only do it when necessary (depending on error types)

    // output cursors or text_positions
    if constexpr (cfg_t::template exists<search_cfg::output<detail::search_output_index_cursor>>())
    {
        return internal_hits;
    }
    else
    {
        std::vector<hit_t> hits;

        if constexpr (cfg_t::template exists<search_cfg::mode<detail::search_mode_best>>())
        {
            // only one cursor is reported but it might contain more than one text position
            if (!internal_hits.empty())
            {
                auto text_pos = internal_hits[0].lazy_locate();
                hits.push_back(text_pos[0]);
            }
            else if (!internal_itv_hits.empty())
            {
                hits.push_back(internal_itv_hits[0]);
            }
        }
        else
        {
            for (auto const & cur : internal_hits)
            {
                for (auto const & text_pos : cur.locate())
                    hits.push_back(text_pos);

//                 std::sort(hits.begin(), hits.end());
//                 hits.erase(std::unique(hits.begin(), hits.end()), hits.end());
            }
            //for more than 4 cursors it is faster to sort once.
            //https://solarianprogrammer.com/2012/10/24/cpp-11-sort-benchmark/
            hits.insert(hits.end(), internal_itv_hits.begin(), internal_itv_hits.end());
            std::sort(hits.begin(), hits.end());
            hits.erase(std::unique(hits.begin(), hits.end()), hits.end());
        }
        return hits;
    }
}

/*!\brief Search a query or a range of queries in an index.
 * \tparam index_t    Must model seqan3::FmIndex.
 * \tparam queries_t  Must model std::ranges::RandomAccessRange over the index's alphabet.
 *                    a range of queries must additionally model std::ranges::ForwardRange.
 * \param[in] index   String index to be searched.
 * \param[in] text   String to be searched. (Optional: Can be empty vector or collection.)
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
template <typename index_t, typename text_t, typename queries_t, typename configuration_t>
inline auto search_all(index_t const & index, text_t const & text, queries_t & queries, configuration_t const & cfg)
{
    using cfg_t = remove_cvref_t<configuration_t>;
    // return type: for each query: a vector of text_positions (or cursors)
    // delegate params: text_position (or cursor). we will withhold all hits of one query anyway to filter
    //                  duplicates. more efficient to call delegate once with one vector instead of calling
    //                  delegate for each hit separately at once.
    using text_pos_t = std::conditional_t<index_t::is_collection_,
                                          std::pair<typename index_t::size_type, typename index_t::size_type>,
                                          typename index_t::size_type>;
    using hit_t = std::conditional_t<cfg_t::template exists<search_cfg::output<detail::search_output_index_cursor>>(),
                                     typename index_t::cursor_type,
                                     text_pos_t>;
    //TODO check if text is collection or not depending on is_collection_
    if constexpr (std::ranges::ForwardRange<queries_t> && std::ranges::RandomAccessRange<value_type_t<queries_t>>)
    {
        // TODO: if constexpr (contains<search_cfg::id::on_hit>(cfg))
        std::vector<std::vector<hit_t>> hits;
        hits.reserve(std::distance(queries.begin(), queries.end()));
        for (auto const query : queries)
        {
            hits.push_back(search_single(index, text, query, cfg));
        }
        return hits;
    }
    else // std::ranges::RandomAccessRange<queries_t>
    {
        // TODO: if constexpr (contains<search_cfg::id::on_hit>(cfg))
        return search_single(index, text, queries, cfg);
    }
}

//!\}

} // namespace seqan3::detail

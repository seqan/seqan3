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

#include <seqan3/core/algorithm/configuration.hpp>
#include <seqan3/range/view/persist.hpp>
#include <seqan3/search/algorithm/detail/search.hpp>
#include <seqan3/search/fm_index/all.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/std/ranges>

namespace seqan3
{

/*!\addtogroup submodule_search_algorithm
 * \{
 */

/*!\brief Search a query or a range of queries in an index.
 * \tparam index_t    Must model seqan3::FmIndex.
 * \tparam queries_t  Must model std::ranges::RandomAccessRange over the index's alphabet.
 *                    a range of queries must additionally model std::ranges::ForwardRange.
 * \param[in] queries A single query or a range of queries.
 * \param[in] index   String index to be searched.
 * \param[in] cfg     A configuration object specifying the search parameters (e.g. number of errors, error types,
 *                    output format, etc.).
 * \returns An object modelling std::ranges::Range containing the hits (the type depends on the specification
            in `cfg`), or `void` if an on_hit delegate has been specified.
 *
 * ### Complexity
 *
 * Each query with \f$e\f$ errors takes \f$O(|query|^e)\f$ where \f$e\f$ is the maximum number of errors.
 *
 * ### Exceptions
 *
 * Strong exception guarantee if iterating the query does not change its state and if invoking a possible delegate
 * specified in `cfg` also has a strong exception guarantee; basic exception guarantee otherwise.
 */
template <FmIndex index_t, typename queries_t, typename configuration_t>
//!\cond
    requires
        (std::ranges::RandomAccessRange<queries_t> ||
            (std::ranges::ForwardRange<queries_t> && std::ranges::RandomAccessRange<value_type_t<queries_t>>)) &&
        detail::is_type_specialisation_of_v<remove_cvref_t<configuration_t>, configuration>
//!\endcond
inline auto search(queries_t && queries, index_t const & index, configuration_t const & cfg)
{
    assert(alphabet_size<innermost_value_type_t<queries_t>> == index.sigma);

    using cfg_t = remove_cvref_t<configuration_t>;

    if constexpr (cfg_t::template exists<search_cfg::max_error>())
    {
        auto & [total, subs, ins, del] = get<search_cfg::max_error>(cfg).value;
        if (subs > total)
            throw std::invalid_argument("The substitution error threshold is higher than the total error threshold.");
        if (ins > total)
            throw std::invalid_argument("The insertion error threshold is higher than the total error threshold.");
        if (del > total)
            throw std::invalid_argument("The deletion error threshold is higher than the total error threshold.");
    }
    else if constexpr (cfg_t::template exists<search_cfg::max_error_rate>())
    {
        auto & [total, subs, ins, del] = get<search_cfg::max_error_rate>(cfg).value;
        if (subs > total)
            throw std::invalid_argument("The substitution error threshold is higher than the total error threshold.");
        if (ins > total)
            throw std::invalid_argument("The insertion error threshold is higher than the total error threshold.");
        if (del > total)
            throw std::invalid_argument("The deletion error threshold is higher than the total error threshold.");
    }

    if constexpr (cfg_t::template exists<search_cfg::mode>())
    {
        if constexpr (cfg_t::template exists<search_cfg::output>())
            return detail::search_all(index, queries, cfg);
        else
            return detail::search_all(index, queries, cfg | search_cfg::output{search_cfg::text_position});
    }
    else
    {
        configuration const cfg2 = cfg | search_cfg::mode{search_cfg::all};
        if constexpr (cfg_t::template exists<search_cfg::output>())
            return detail::search_all(index, queries, cfg2);
        else
            return detail::search_all(index, queries, cfg2 | search_cfg::output{search_cfg::text_position});
    }
}

//! \overload
template <FmIndex index_t, typename configuration_t>
inline auto search(char const * const queries, index_t const & index, configuration_t const & cfg)
{
    return search(std::string_view{queries}, index, cfg);
}

//! \overload
template <FmIndex index_t, typename configuration_t>
inline auto search(std::initializer_list<char const * const> const & queries,
                   index_t const & index,
                   configuration_t const & cfg)
{
    std::vector<std::string_view> query;
    query.reserve(std::ranges::size(queries));
    std::ranges::for_each(queries, [&query] (char const * const q) { query.push_back(std::string_view{q}); });
    return search(query, index, cfg);
}

/*!\brief Search a query or a range of queries in an index.
 *        It will not allow for any errors and will output all matches as positions in the text.
 * \tparam index_t    Must model seqan3::FmIndex.
 * \tparam queries_t  Must model std::ranges::RandomAccessRange over the index's alphabet.
 *                    a range of queries must additionally model std::ranges::ForwardRange.
 * \param[in] queries A single query or a range of queries.
 * \param[in] index   String index to be searched.
 * \returns An object modelling std::ranges::Range containing the hits as positions in the searched text.
 *
 * ### Complexity
 *
 * Each query with \f$e\f$ errors takes \f$O(|query|^e)\f$ where \f$e\f$ is the maximum number of errors.
 *
 * ### Exceptions
 *
 * Strong exception guarantee if iterating the query does not change its state; basic exception guarantee otherwise.
 */
template <FmIndex index_t, typename queries_t>
//!\cond
    requires std::ranges::RandomAccessRange<queries_t> ||
             (std::ranges::ForwardRange<queries_t> && std::ranges::RandomAccessRange<value_type_t<queries_t>>)
//!\endcond
inline auto search(queries_t && queries, index_t const & index)
{
    assert(alphabet_size<innermost_value_type_t<queries_t>> == index.sigma);

    configuration const default_cfg = search_cfg::max_error{search_cfg::total{0},
                                                            search_cfg::substitution{0},
                                                            search_cfg::insertion{0},
                                                            search_cfg::deletion{0}}
                                            | search_cfg::output{search_cfg::text_position}
                                            | search_cfg::mode{search_cfg::all};
    return search(queries, index, default_cfg);
}

//! \overload
template <FmIndex index_t>
inline auto search(index_t const & index, char const * const queries)
{
    return search(std::string_view{queries}, index);
}

//! \overload
template <FmIndex index_t>
inline auto search(std::initializer_list<char const * const> const & queries, index_t const & index)
{
    std::vector<std::string_view> query;
    query.reserve(std::ranges::size(queries));
    std::ranges::for_each(queries, [&query] (char const * const q) { query.push_back(std::string_view{q}); });
    return search(query, index);
}

//!\}

} // namespace seqan3

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
#include <seqan3/range/views/persist.hpp>
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
 * \tparam index_t    Must model seqan3::fm_index_specialisation.
 * \tparam queries_t  Must model std::ranges::random_access_range over the index's alphabet and std::ranges::sized_range.
 *                    A range of queries must additionally model std::ranges::forward_range and std::ranges::sized_range.
 * \param[in] queries A single query or a range of queries.
 * \param[in] index   String index to be searched.
 * \param[in] cfg     A configuration object specifying the search parameters (e.g. number of errors, error types,
 *                    output format, etc.).
 * \returns
 * <table>
 *   <tr>
 *     <th>seqan3::text_layout</th>
 *     <th>seqan3::search_cfg::output</th>
 *     <th>Result of seqan3::search()</th>
 *   </tr>
 *   <tr>
 *     <td style="text-align:center">\ref seqan3::text_layout "single"</td>
 *     <td style="text-align:center">\ref seqan3::search_cfg::text_position "text_position"</td>
 *     <td>A `std::vector<size_t>` representing text positions where the search was successful.</td>
 *   </tr>
 *   <tr>
 *     <td style="text-align:center">\ref seqan3::text_layout "single"</td>
 *     <td style="text-align:center">\ref seqan3::search_cfg::index_cursor "index_cursor"</td>
 *     <td>A `std::vector<typename index_t::cursor_type>` containing index_cursors at the text positions where the
 *         search was successful.</td>
 *   </tr>
 *   <tr>
 *     <td style="text-align:center">\ref seqan3::text_layout "collection"</td>
 *     <td style="text-align:center">\ref seqan3::search_cfg::text_position "text_position"</td>
 *     <td>A `std::vector<std::pair<size_t, size_t>>` where the first element of the
 *         `std::pair` specifies the text index in the collection and the second element contains the position in that
 *         text of the collection.</td>
 *   </tr>
 *   <tr>
 *     <td style="text-align:center">\ref seqan3::text_layout "collection"</td>
 *     <td style="text-align:center">\ref seqan3::search_cfg::index_cursor "index_cursor"</td>
 *     <td>A `std::vector<typename index_t::cursor_type>` containing index_cursors at the text positions where the
 *         search was successful.</td>
 *   </tr>
 * </table>
 *
 * \if DEV \note Always returns `void` if an on_hit delegate has been specified.\endif
 *
 * \details
 *
 * **Header**
 * ```cpp
 *      #include <seqan3/search/algorithm/search.hpp>
 * ```
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
template <fm_index_specialisation index_t, typename queries_t, typename configuration_t = decltype(search_cfg::default_configuration)>
inline auto search(queries_t && queries,
                   index_t const & index,
                   configuration_t const & cfg = search_cfg::default_configuration)
{
    if constexpr(dimension_v<queries_t> == 1u)
    {
        static_assert(std::ranges::random_access_range<queries_t>, "The query sequence must model random_access_range.");
        static_assert(std::ranges::sized_range<queries_t>, "The query sequence must model sized_range.");
    }
    else
    {
        static_assert(std::ranges::forward_range<queries_t>, "The query collection must model forward_range.");
        static_assert(std::ranges::sized_range<queries_t>, "The query collection must model sized_range.");
        static_assert(std::ranges::random_access_range<value_type_t<queries_t>>,
                      "Elements of the query collection must model random_access_range.");
        static_assert(std::ranges::sized_range<value_type_t<queries_t>>,
                      "Elements of the query collection must model sized_range.");
    }

    static_assert(detail::is_type_specialisation_of_v<remove_cvref_t<configuration_t>, configuration>,
                  "cfg must be a specialisation of seqan3::configuration.");

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
//!\cond DEV
//! \overload
template <fm_index_specialisation index_t, typename configuration_t = decltype(search_cfg::default_configuration)>
inline auto search(char const * const queries,
                   index_t const & index,
                   configuration_t const & cfg = search_cfg::default_configuration)
{
    return search(std::string_view{queries}, index, cfg);
}

//! \overload
template <fm_index_specialisation index_t, typename configuration_t = decltype(search_cfg::default_configuration)>
inline auto search(std::initializer_list<char const * const> const & queries,
                   index_t const & index,
                   configuration_t const & cfg = search_cfg::default_configuration)
{
    std::vector<std::string_view> query;
    query.reserve(std::ranges::size(queries));
    std::ranges::for_each(queries, [&query] (char const * const q) { query.push_back(std::string_view{q}); });
    return search(query, index, cfg);
}
//!\endcond

//!\}

} // namespace seqan3

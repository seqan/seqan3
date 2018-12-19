// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ============================================================================

/*!\file
 * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
 * \brief Provides the public interface for search algorithms.
 */

#pragma once

#include <seqan3/core/algorithm/configuration.hpp>
#include <seqan3/range/view/persist.hpp>
#include <seqan3/search/algorithm/detail/search.hpp>
#include <seqan3/search/fm_index/all.hpp>

namespace seqan3
{

/*!\addtogroup submodule_search_algorithm
 * \{
 */

/*!\brief Search a query or a range of queries in an index.
 * \tparam index_t    Must model seqan3::FmIndex.
 * \tparam queries_t  Must be a std::ranges::RandomAccessRange over the index's alphabet.
 *                    a range of queries must additionally model std::ranges::ForwardRange.
 * \param[in] index   String index to be searched.
 * \param[in] queries A single query or a range of queries.
 * \param[in] cfg     A configuration object specifying the search parameters (e.g. number of errors, error types,
 *                    output format, etc.).
 * \returns An object modelling std::ranges::Range containing the hits (the type depends on the specification
            in `cfg`), or `void` if an on_hit delegate has been specified.
 *
 * \todo Update concepts and documentation of `configuration_t` everywhere once it has been refactored by rrahn.
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
inline auto search(index_t const & index, queries_t && queries, configuration_t const & cfg)
{
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

/*!\brief Search a query or a range of queries in an index.
 *        It will not allow for any errors and will output all matches as positions in the text.
 * \tparam index_t    Must model seqan3::FmIndex.
 * \tparam queries_t  Must be a std::ranges::RandomAccessRange over the index's alphabet.
 *                    a range of queries must additionally model std::ranges::ForwardRange.
 * \param[in] index   String index to be searched.
 * \param[in] queries A single query or a range of queries.
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
inline auto search(index_t const & index, queries_t && queries)
{
    configuration const default_cfg = search_cfg::max_error{search_cfg::total{0},
                                                            search_cfg::substitution{0},
                                                            search_cfg::insertion{0},
                                                            search_cfg::deletion{0}}
                                            | search_cfg::output{search_cfg::text_position}
                                            | search_cfg::mode{search_cfg::all};
    return search(index, queries, default_cfg);
}

//!\}

} // namespace seqan3

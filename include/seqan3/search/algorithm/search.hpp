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
 * \brief
 */

#pragma once

#include <seqan3/range/view/persist.hpp>
#include <seqan3/search/algorithm/detail/search.hpp>
#include <seqan3/search/fm_index/all.hpp>

namespace seqan3
{

//!\brief \todo Document!
template <typename index_t, typename queries_t, typename config_t>
//!\cond
    requires
        (std::ranges::RandomAccessRange<queries_t> ||
            (std::ranges::ForwardRange<queries_t> && std::ranges::RandomAccessRange<value_type_t<queries_t>>)) &&
        detail::is_algorithm_configuration_v<remove_cvref_t<config_t>>
//!\endcond
inline auto search(index_t const & index, queries_t const & queries, config_t const & cfg)
{
    // TODO: replace enumeration of all code paths to set all required configuration objects
    if constexpr (contains<search_cfg::id::mode>(cfg))
    {
        if constexpr (contains<search_cfg::id::output>(cfg))
            return detail::_search(index, queries, cfg);
        else
            return detail::_search(index, queries, cfg | search_cfg::output(search_cfg::text_position));
    }
    else
    {
        // TODO: overload pipe operator for empty config object
        if constexpr (std::Same<remove_cvref_t<decltype(cfg)>, detail::configuration<>>)
        {
            detail::configuration const cfg2 = search_cfg::mode(search_cfg::all);
            if constexpr (contains<search_cfg::id::output>(cfg))
                return detail::_search(index, queries, cfg2);
            else
                return detail::_search(index, queries, cfg2 | search_cfg::output(search_cfg::text_position));
        }
        else
        {
            detail::configuration const cfg2 = cfg | search_cfg::mode(search_cfg::all);
            if constexpr (contains<search_cfg::id::output>(cfg))
                return detail::_search(index, queries, cfg2);
            else
                return detail::_search(index, queries, cfg2 | search_cfg::output(search_cfg::text_position));
        }
    }
}

// TODO: const &, &&, etc.: use forwarding references and add a view::persist
// DOC: insertion/deletion are with resp. to the query. i.e. an insertion is the insertion of a base into the query
// that does not occur in the text at the position
//!\brief \todo Document!
template <typename index_t, typename queries_t>
//!\cond
    requires std::ranges::RandomAccessRange<queries_t> ||
             (std::ranges::ForwardRange<queries_t> && std::ranges::RandomAccessRange<value_type_t<queries_t>>)
//!\endcond
inline auto search(index_t const & index, queries_t const & queries)
{
    // TODO: auto const queries_lvalue = queries | view::persist;

    detail::configuration const default_cfg = search_cfg::max_error(search_cfg::total{0}, search_cfg::substitution{0},
                                                            search_cfg::insertion{0}, search_cfg::deletion{0})
                                            | search_cfg::output(search_cfg::text_position)
                                            | search_cfg::mode(search_cfg::all);
    return search(index, queries/*_lvalue*/, default_cfg);
}

} // namespace seqan3

// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 * \brief Provides the seqan3::detail::policy_max_error.
 */

#pragma once

#include <seqan3/std/ranges>

#include <seqan3/search/detail/search_common.hpp>
#include <seqan3/search/detail/search_traits.hpp>

namespace seqan3::detail
{

//!\brief Provides the function `max_error_counts` if inherited by a search algorithm.
//!\ingroup search
struct policy_max_error
{
protected:
    /*!\brief Returns a detail::search_param object filled by the information from the configuration.
     * \tparam configuration_t The search configuration type.
     * \tparam query_t Must model std::ranges::foward_range over the index's alphabet.
     * \param[in] cfg The configuration object.
     * \param[in] query The current query sequence.
     *
     * \details
     * 
     * If seqan3::max_error is set in the configuration, its value already is of type detail::search_param and
     * can be returned as is.
     * If seqan3::max_error_rate is set in the configuration, the error rates are converted to error counts
     * based on the length of the query sequence.
     */
    template <typename configuration_t, std::ranges::forward_range query_t>
    auto max_error_counts(configuration_t const & cfg, [[maybe_unused]] query_t && query)
    {
        using search_traits_t = search_traits<configuration_t>;

        // retrieve error numbers / rates
        detail::search_param max_error{0, 0, 0, 0};

        if constexpr (search_traits_t::search_with_max_error)
        {
            return get<search_cfg::max_error>(cfg).value;
        }
        else if constexpr (search_traits_t::search_with_max_error_rate)
        {
            // NOTE: Casting doubles rounds towards zero (i.e. floor for positive numbers). Thus, given a rate of
            // 10% and a read length of 101 the max number of errors is correctly casted from 10.1 errors to 10
            auto & [total, subs, ins, del] = max_error;
            std::tie(total, subs, ins, del) = std::apply([&query] (auto && ... args)
            {
                return std::tuple{(args * std::ranges::size(query))...};
            }, get<search_cfg::max_error_rate>(cfg).value);
        }

        return max_error;
    }
};

} // namespace seqan3::detail

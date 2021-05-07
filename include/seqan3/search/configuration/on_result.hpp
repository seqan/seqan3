// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::search_cfg::on_result.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <type_traits>


#include <seqan3/core/algorithm/pipeable_config_element.hpp>
#include <seqan3/core/semiregular_box.hpp>
#include <seqan3/search/configuration/detail.hpp>

namespace seqan3::search_cfg
{

/*!\brief Configuration element to provide a user defined callback function for the search.
 * \ingroup search_configuration
 *
 * \tparam callback_t The type of the callback; must model std::invocable with the generated seqan3::search_result
 *                    and std::move_constructible.
 *
 * \details
 *
 * Allows the user to specify a callback that should be called for every computed search result. The callback
 * must take exactly one argument for the search result and return `void`. If the user callback is
 * specified, the call to the search algorithm seqan3::search will return nothing, i.e. it does not return
 * a seqan3::algorithm_result_generator_range any more. Note that within a parallel configuration, the order of the
 * generated search results and therefore the call to the user callback is non-deterministic.
 * However, the continuation interface with the
 * user callback can be more efficient in a concurrent environment.
 *
 * \if DEV
 * The given callback is wrapped inside a seqan3::semiregular_box wrapper type. This allows to also
 * use lambdas with a capture block, which otherwise are not std::copy_assignable and therefore invalidate the
 * requirements for the configuration element (must model std::semiregular).
 * \endif
 *
 * ### Example
 *
 * The following code snippet demonstrates the basic usage:
 *
 * \include test/snippet/search/configuration_on_result.cpp
 */
template <std::move_constructible callback_t>
class on_result : public seqan3::pipeable_config_element
{
public:
    //!\brief The stored callable which will be invoked with the search result.
    seqan3::semiregular_box_t<callback_t> callback{}; // Allow lambdas with capture block which are not copy_assignable.

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr on_result() = default; //!< Defaulted.
    constexpr on_result(on_result const &) = default; //!< Defaulted.
    constexpr on_result(on_result &&) = default; //!< Defaulted.
    constexpr on_result & operator=(on_result const &) = default; //!< Defaulted.
    constexpr on_result & operator=(on_result &&) = default; //!< Defaulted.
    ~on_result() = default; //!< Defaulted.

    /*!\brief Constructs the configuration element with the given user callback.
     * \param[in] callback The callback to invoke with a computed seqan3::search_result.
     */
    constexpr explicit on_result(callback_t callback) : callback{std::forward<callback_t>(callback)}
    {}
    //!\}

    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr seqan3::detail::search_config_id id{seqan3::detail::search_config_id::on_result};
};

/*!\name Type deduction guides
 * \{
 */
//!\brief Deduces the callback type from a forwarding constructor argument.
template <std::move_constructible callback_t>
on_result(callback_t &&) -> on_result<std::decay_t<callback_t>>;
//!\}
}  // namespace seqan3::search_cfg

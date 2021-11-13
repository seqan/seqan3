// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 * \brief Provides the concept seqan3::detail::is_container_option.
 */

#pragma once

#include <seqan3/std/concepts>
#include <string>
#include <seqan3/std/type_traits>

#include <seqan3/core/platform.hpp>

namespace seqan3::detail
{

/*!\interface seqan3::detail::is_container_option <>
 * \ingroup argument_parser
 * \brief Decides whether an option type is considered a list.
 * \details
 *
 * In the seqan3::argument_parser::add_option or seqan3::argument_parser::add_positional_argument calls,
 * it is important whether the type of the value, in which the respective option is stored (option_type),
 * is a container or not because nternal handling of option containers is different.
 *
 * This concept decides whether an option type is considered a container or not.
 * In general, all standard library containers except std::string can be considered option containers.
 *
 * In order to be considered an option container, the `option_type` as to fulfil the following:
 * * `option_type` may not be `std::string`
 * * Member type `option_type::value_type` must exist
 * * Member function `option_type::push_back(typename option_type::value_type) -> void` must exist.
 *
 * \noapi{Exposition only.}
 */
//!\cond
template <typename option_type>
concept is_container_option = !std::is_same_v<std::remove_cvref_t<option_type>, std::string> &&
                                     requires (option_type container,
                                               typename std::remove_reference_t<option_type>::value_type value)
{
    { container.push_back(value) };
};
//!\endcond

} // namespace seqan3::detail

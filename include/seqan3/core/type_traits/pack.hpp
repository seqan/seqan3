// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides unary type traits on a set of types, usually provided as template argument pack.
 */

#pragma once

#include <type_traits>

#include <seqan3/core/platform.hpp>

namespace seqan3::detail
{

//!\addtogroup type_traits
//!\{

//!\brief Indicates whether the first template argument is contained in the remaining.
//!\implements seqan3::unary_type_trait
template <typename target_t, typename ...pack>
struct type_in_pack : std::false_type {};

//!\cond
template <typename target_t, typename ...pack>
struct type_in_pack<target_t, target_t, pack...> : std::true_type {};

template <typename target_t, typename pack1, typename ...pack>
struct type_in_pack<target_t, pack1, pack...> : type_in_pack<target_t, pack...> {};
//!\endcond

//!\brief Shortcut for seqan3::detail::type_in_pack (unary_type_trait shortcut).
//!\relates seqan3::detail::type_in_pack
template <typename target_t, typename ...pack>
inline bool constexpr type_in_pack_v = type_in_pack<target_t, pack...>::value;

//!\}

} // namespace seqan3::detail

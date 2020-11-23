// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief [DEPRECATED] Provides unary type traits on a set of types, usually provided as template argument pack.
 * \deprecated This header is deprecated and will be removed in SeqAn-3.1.
 *             Please use seqan3::pack_traits::contains in <seqan3/core/type_list/traits.hpp> instead.
 */

#pragma once

#include <seqan3/core/platform.hpp>
#include <seqan3/core/type_list/traits.hpp>

SEQAN3_DEPRECATED_HEADER("This header is deprecated and will be removed in SeqAn-3.1. Please use seqan3::pack_traits::contains in <seqan3/core/type_list/traits.hpp> instead.")

namespace seqan3::detail
{

//!\addtogroup type_traits
//!\{

/*!\brief [DEPRECATED] Indicates whether the first template argument is contained in the remaining.
 * \implements seqan3::unary_type_trait
 * \deprecated This struct is deprecated and will be removed in SeqAn-3.1.
 *             Please use seqan3::pack_traits::contains instead.
 */
template <typename target_t, typename ...pack>
struct SEQAN3_DEPRECATED_310 type_in_pack : std::conditional_t<seqan3::pack_traits::contains<target_t, pack...>,
                                                               std::true_type,
                                                               std::false_type> {};

/*!\brief [DEPRECATED] Shortcut for seqan3::detail::type_in_pack (unary_type_trait shortcut).
 * \deprecated This variable is deprecated and will be removed in SeqAn-3.1.
 *             Please use seqan3::pack_traits::contains instead.
 */
template <typename target_t, typename ...pack>
SEQAN3_DEPRECATED_310 inline bool constexpr type_in_pack_v = seqan3::pack_traits::contains<target_t, pack...>;

//!\}

} // namespace seqan3::detail

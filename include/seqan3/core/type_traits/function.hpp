// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides various type traits for use on functions.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/platform.hpp>

// ----------------------------------------------------------------------------
// is_constexpr
// ----------------------------------------------------------------------------

/*!\brief Returns true if the expression passed to this macro can be evaluated at compile time, false otherwise.
 * \ingroup type_traits
 * \returns true or false.
 */
#define SEQAN3_IS_CONSTEXPR(...) std::integral_constant<bool, __builtin_constant_p((__VA_ARGS__, 0))>::value

// ----------------------------------------------------------------------------
// multi_invocable
// ----------------------------------------------------------------------------

namespace seqan3::detail
{

/*!\brief A type that can conveniently inherit multiple invocables and acts as a union over them.
 * \tparam invocable_ts The types to inherit from.
 * \ingroup type_traits
 */
template <typename ...invocable_ts>
struct multi_invocable : invocable_ts...
{
    //!\brief Inherit the function call operators.
    using invocable_ts::operator()...;
};

//!\brief Deduction guides for seqan3::detail::multi_invocable.
template <typename ...invocable_ts>
multi_invocable(invocable_ts...) -> multi_invocable<invocable_ts...>;

} // namespace seqan3::detail

// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides concept for seqan3::function_like concept.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <functional>
#include <type_traits>

#include <seqan3/core/platform.hpp>

namespace seqan3
{

// ----------------------------------------------------------------------------
// function_like
// ----------------------------------------------------------------------------

/*!\interface seqan3::function_like
 * \ingroup core
 * \brief Checks if a type is a function type, a lambda function type, a std::function type or a function pointer type.
 *
 * \tparam t The type to check if it is a valid function, lambda function, std::function or function pointer type.
 *
 * \details
 *
 * This concept checks if the given type is a lambda function type, a std::function type or a function pointer type.
 *
 * ### Example
 *
 * \include test/snippet/core/concept/function_like.cpp
 */
/*!\name Requirements for seqan3::function_like
 * \brief You can expect these functions on all types that implement seqan3::function_like.
 * \{
 */
/*!\fn              auto && std::get<i>(type && val)
 * \brief           Return the i-th element of the tuple.
 * \relates         seqan3::function_like
 * \tparam          i The index of the element to return (of type `size_t`).
 * \param[in,out]   val The tuple-like object to operate on.
 * \returns         The i-th value in the tuple.
 *
 * \details
 * \attention This is a concept requirement, not an actual function (however types satisfying this concept
 * will provide an implementation).
 * \attention This constraint is **not enforced** since empty tuples are valid.
 */
template <typename t>
SEQAN3_CONCEPT function_like = std::is_function_v<t> || requires (t v ) { {std::function{v}}; };

} // namespace seqan3

// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin 
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik 
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides various metafunctions on generic types.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <type_traits>

#include <range/v3/utility/iterator_traits.hpp>
#include <range/v3/range_traits.hpp>

#include <seqan3/core/platform.hpp>

namespace seqan3
{

/*!\addtogroup metafunction
 * \{
 */

// ----------------------------------------------------------------------------
// remove_cvref_t
// ----------------------------------------------------------------------------

/*!\brief Return the input type with `const`, `volatile` and references removed [Type metafunction].
 * \tparam t The type to operate on.
 */
template <typename t>
using remove_cvref_t = std::remove_cv_t<std::remove_reference_t<t>>;

//!\}

} // namespace seqan3

namespace seqan3::detail
{

/*!\addtogroup metafunction
 * \{
 */

// ----------------------------------------------------------------------------
// remove_cvref_t
// ----------------------------------------------------------------------------

//!\brief Return the type of std::ignore with `const`, `volatile` and references removed [Type metafunction].
using ignore_t = remove_cvref_t<decltype(std::ignore)>;

/*!\brief Return whether the input type with `const`, `volatile` and references removed is std::ignore's type.
 * [Value metafunction].
 * \tparam t The type to operate on.
 */
template <typename t>
constexpr bool decays_to_ignore_v = std::is_same_v<remove_cvref_t<t>, ignore_t>;

//!\}

} // namespace seqan3::detail

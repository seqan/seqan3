// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::semiregular_box.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/detail/copyable_wrapper.hpp>

SEQAN3_DEPRECATED_HEADER("This header is deprecated and will be removed in SeqAn-3.3.0; Please #include "
                         "<seqan3/core/detail/copyable_wrapper.hpp> instead.")

namespace seqan3
{
/*!\brief Utility wrapper that behaves like std::optional but makes the type conform with the std::semiregular concept.
 * \ingroup core
 *
 * \see https://en.cppreference.com/w/cpp/ranges/copyable_wrapper
 *
 * \deprecated Please use seqan3::detail::copyable_wrapper
 */
template <typename t>
using semiregular_box SEQAN3_DEPRECATED_330 = detail::copyable_wrapper<t>;

/*!\brief Utility transformation trait to get a wrapper type that models std::semiregular.
 * \ingroup core
 *
 * \see https://en.cppreference.com/w/cpp/ranges/copyable_wrapper
 *
 * \deprecated Please use seqan3::detail::copyable_wrapper
 */
template <typename t>
using semiregular_box_t SEQAN3_DEPRECATED_330 = detail::copyable_wrapper_t<t>;

} // namespace seqan3

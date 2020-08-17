// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::semiregular_box.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <range/v3/utility/semiregular_box.hpp>

#include <seqan3/core/platform.hpp>

namespace seqan3
{
/*!\brief Utility wrapper that behaves like std::optional but makes the type conform with the std::semiregular concept.
          Imported from ranges::semiregular_box.
 * \ingroup core
 *
 * \see https://en.cppreference.com/w/cpp/ranges/semiregular_wrapper
 */
using SEQAN3_DOXYGEN_ONLY(semiregular_box =) ::ranges::semiregular_box;

/*!\brief Utility transformation trait to get a wrapper type that models std::semiregular. Imported from
          ranges::semiregular_box_t.
 * \ingroup core
 *
 * \see https://en.cppreference.com/w/cpp/ranges/semiregular_wrapper
 */
using SEQAN3_DOXYGEN_ONLY(semiregular_box_t =) ::ranges::semiregular_box_t;

}  // namespace seqan3

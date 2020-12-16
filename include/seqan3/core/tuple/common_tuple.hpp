// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::common_tuple and seqan3::common_pair.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <range/v3/utility/common_tuple.hpp>

#include <seqan3/core/platform.hpp>

namespace seqan3
{

/*!\brief A common tuple type that behaves like a regular std::tuple, but can be used as a reference type proxy for
 *        output iterators.
 *
 * \details
 *
 * Alias definition of the ranges::common_tuple.
 */
using SEQAN3_DOXYGEN_ONLY(common_tuple =) ::ranges::common_tuple;

/*!\brief A common pair type that behaves like a regular std::pair, but can be used as a reference type proxy for
 *        output iterators.
 *
 * \details
 *
 * Alias definition of the ranges::common_pair.
 */
using SEQAN3_DOXYGEN_ONLY(common_pair =) ::ranges::common_pair;

}  // namespace seqan3

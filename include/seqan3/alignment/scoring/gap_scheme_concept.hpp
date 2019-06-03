// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::GapScheme.
 * \author René Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <type_traits>

#include <seqan3/std/concepts>

namespace seqan3
{

/*!\interface seqan3::GapScheme <>
 * \brief A concept that requires that the type gives access to gap penalties.
 * \tparam t            The type the concept check is performed on (the putative gap scheme).
 * \ingroup scoring
 *
 * \details
 *
 * This concept makes no assumptions about configurability or assignability of the gap scheme, only the
 * ability to request the gap penalties.
 *
 */

/*!\name Requirements for seqan3::GapScheme
 * \brief You can expect these members on all types that implement seqan3::GapScheme.
 * \memberof seqan3::GapScheme
 * \{
 */

/*!\fn       ptrdiff_t score(size_t s);
 * \brief    Returns the score for a gap of length `s`.
 *
 * \details
 * \attention This is a concept requirement, not an actual function (however types satisfying this concept
 * will provide an implementation).
 */
//!\}

//!\cond
template <typename t>
SEQAN3_CONCEPT GapScheme = requires (t scheme, size_t num)
{
    { scheme.score(num) } -> ptrdiff_t;
};
//!\endcond

} // namespace seqan3

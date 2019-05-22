// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::ostream and seqan3::ostreambuf iterator.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <iterator>

#ifndef __cpp_lib_ranges
#include <range/v3/iterator/stream_iterators.hpp>
#endif // __cpp_lib_ranges

#include <seqan3/core/platform.hpp>

namespace seqan3
{

/*!\addtogroup stream
 * \{
 */

#ifdef __cpp_lib_ranges

/*!\typedef seqan3::ostream_iterator
* \brief Alias for std::ostream_iterator. Writes successive elements onto the output stream from which it was
*        constructed.
*/
using SEQAN3_DOXYGEN_ONLY(ostream_iterator =) ::std::ostream_iterator;

/*!\typedef seqan3::ostreambuf_iterator
 * \brief Alias for std::ostreambuf_iterator. Writes successive characters onto the output stream from which it
 *        was constructed.
 */
using SEQAN3_DOXYGEN_ONLY(ostreambuf_iterator =) ::std::ostreambuf_iterator;

#else

/*!\typedef seqan3::ostream_iterator
 * \brief Alias for ranges::ostream_iterator. Writes successive elements onto the output stream from which it was
 *        constructed.
 */
using SEQAN3_DOXYGEN_ONLY(ostream_iterator =) ::ranges::ostream_iterator;

/*!\typedef seqan3::ostreambuf_iterator
 * \brief Alias for ranges::ostreambuf_iterator. Writes successive characters onto the output stream from which it
 *        was constructed.
 */
using SEQAN3_DOXYGEN_ONLY(ostreambuf_iterator =) ::ranges::ostreambuf_iterator;

#endif // __cpp_lib_ranges
//!\}
} // namespace seqan3

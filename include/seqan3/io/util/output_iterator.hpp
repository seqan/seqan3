// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2017, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2017, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ============================================================================

/*!\file
 * \ingroup io_util
 * \brief Convenience functions to create single pass output iterator for various source types.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <type_traits>

#include <seqan3/io/concept.hpp>
#include <seqan3/range/container/concept.hpp>

#include <range/v3/utility/associated_types.hpp>
#include <range/v3/utility/iterator.hpp>

namespace seqan3
{
/*!\brief A convenience function template that constructs a ranges::back_insert_iterator for the container \c c with
 * the type deduced from the type of the argument.
 * \tparam container_type Any type that satisfies the seqan3::sequence_concept.
 * \param c The instance to construct the back_insert_iterator for.
 * \ingroup io_util
 */
template <typename container_type>
//!\cond
    requires sequence_concept<std::remove_reference_t<container_type>> &&
             !std::is_const_v<std::remove_reference_t<container_type>>
//!\endcond
auto single_pass_output_iterator(container_type & c)
{
    return ranges::back_insert_iterator<container_type>{c};
}

/*!\brief A convenience function template that constructs an ranges::ostreambuf_iterator for the ostream \c s with
 * the char type and traits type deduced from the type of the argument.
 * \tparam ostream_type Any type that satisfies the seqan3::ostream_concept.
 * \param s The stream instance to construct the ostreambuf_iterator for.
 * \ingroup io_util
 */
template <typename ostream_type>
//!\cond
    requires ostream_concept<std::remove_reference_t<ostream_type>,
                             typename ranges::value_type<std::remove_reference_t<ostream_type>>::type>
//!\endcond
auto single_pass_output_iterator(ostream_type & s)
{
    using char_type   = typename std::remove_reference_t<ostream_type>::char_type;
    using traits_type = typename std::remove_reference_t<ostream_type>::traits_type;
    return ranges::ostreambuf_iterator<char_type, traits_type>{s};
}
} // namespace seqan3

// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
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
 * \brief Provides tokenisation functionality.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <functional>

#include <seqan3/core/concept/core.hpp>
#include <seqan3/core/concept/iterator.hpp>
#include <seqan3/core/metafunction/range.hpp>
#include <seqan3/io/exception.hpp>
#include <seqan3/range/concept.hpp>

namespace seqan3::detail
{
// ----------------------------------------------------------------------------
// transfer_data
// ----------------------------------------------------------------------------

//!\cond
//TODO(rrahn): Replace by anonymous lambda expression once https://gcc.gnu.org/bugzilla/show_bug.cgi?id=85513 is fixed.
auto invocable_dummy = [](auto){};
//!\endcond

template <typename            out_iterator_type,
          input_range_concept input_rng_type,
          typename            delimiter_type,
          typename            asserter_type  = decltype(invocable_dummy)>
//!\cond
    requires output_iterator_concept<std::remove_reference_t<out_iterator_type>, char> &&
             predicate_concept<std::remove_reference_t<delimiter_type>, char> &&
             invocable_concept<std::remove_reference_t<asserter_type>, char>
//!\endcond
void transfer_data(out_iterator_type && receiver,
                   input_rng_type    && transmitter,
                   delimiter_type    && delim,
                   asserter_type     && asserter = std::move(invocable_dummy))
{
    for (auto && c : transmitter)
    {
        if (delim(c)) // delimiter was reached.
            return;
        *receiver = (asserter(c), c);  // Check if c satisfies given assert condition and write to the receiver on success.
        ++receiver;
    }
    throw unexpected_end_of_error{"Reached end of input while expecting more data."};
}

} // namespace seqan3::detail

namespace seqan3
{
// std::interface of istream
// get([std::streamsize count][, char_type delim])
// get_line([delim])
// ignore(count, delim)
// read(count)

// std::interface of ostream
// put
// write
} // namespace seqan3

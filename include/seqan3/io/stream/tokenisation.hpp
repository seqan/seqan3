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
#include <tuple>

#include <seqan3/core/concept/core.hpp>
#include <seqan3/core/concept/iterator.hpp>
#include <seqan3/core/metafunction/basic.hpp>
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
          typename            asserter_type  = decltype(invocable_dummy) &>
//!\cond
    requires output_iterator_concept<std::remove_reference_t<out_iterator_type>, char> &&
             predicate_concept<std::remove_reference_t<delimiter_type>, char> &&
             invocable_concept<std::remove_reference_t<asserter_type>, char>
//!\endcond
void transfer_data(out_iterator_type && receiver,
                   input_rng_type    && transmitter,
                   delimiter_type    && delim,
                   asserter_type     && asserter = invocable_dummy)
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

// ----------------------------------------------------------------------------
// read_until
// ----------------------------------------------------------------------------

template <typename receiver_type,
          input_range_concept input_rng_type,
          typename delimiter_type,
          typename asserter_type  = decltype(detail::invocable_dummy) &>
//!\cond
    requires (std::is_same_v<remove_cvref_t<receiver_type>, remove_cvref_t<decltype(std::ignore)>> ||
              output_iterator_concept<std::remove_reference_t<receiver_type>, char>) &&
             predicate_concept<std::remove_reference_t<delimiter_type>, char> &&
             invocable_concept<std::remove_reference_t<asserter_type>, char>
//!\endcond
void read_until(receiver_type  && rcvr,
                input_rng_type && input_rng,
                delimiter_type && delim,
                asserter_type  && asserter = detail::invocable_dummy)
{
    if constexpr (std::is_same_v<remove_cvref_t<receiver_type>, remove_cvref_t<decltype(std::ignore)>>)
        detail::transfer_data(detail::make_conversion_output_iterator(rcvr),
                              std::forward<input_rng_type>(input_rng),
                              std::forward<delimiter_type>(delim),
                              std::forward<asserter_type>(asserter));
    else
        detail::transfer_data(std::forward<receiver_type>(rcvr),
                              std::forward<input_rng_type>(input_rng),
                              std::forward<delimiter_type>(delim),
                              std::forward<asserter_type>(asserter));
}

// ----------------------------------------------------------------------------
// read_line
// ----------------------------------------------------------------------------

template <typename            receiver_type,
          input_range_concept input_rng_type,
          typename            asserter_type  = decltype(detail::invocable_dummy) &>
//!\cond
    requires (std::is_same_v<remove_cvref_t<receiver_type>, remove_cvref_t<decltype(std::ignore)>> ||
              output_iterator_concept<std::remove_reference_t<receiver_type>, char>) &&
              invocable_concept<std::remove_reference_t<asserter_type>, char>
//!\endcond
void read_line(receiver_type  && rcvr,
               input_rng_type && input_rng,
               asserter_type  && asserter = detail::invocable_dummy)
{
    read_until(std::forward<receiver_type>(rcvr),
               std::forward<input_rng_type>(input_rng),
               is_char<'\n'>{} || is_char<'\r'>{},
               std::forward<asserter_type>(asserter));

    // Check if the statement was carriage return plus new-line.
    auto it = ranges::begin(input_rng);
    if (*it == '\r')
        if (*(++it) != '\n')  // consume the '\r' symbol and check if '\n' follows.
            throw parse_error{"Missing newline '\n' character after reading '\r' character."};
    ++it; // extract the newline character.
}

// ----------------------------------------------------------------------------
// read_n
// ----------------------------------------------------------------------------

template <typename            receiver_type,
          input_range_concept input_rng_type,
          typename            asserter_type  = decltype(detail::invocable_dummy) &>
//!\cond
    requires (std::is_same_v<remove_cvref_t<receiver_type>, remove_cvref_t<decltype(std::ignore)>> ||
              output_iterator_concept<std::remove_reference_t<receiver_type>, char>) &&
              invocable_concept<std::remove_reference_t<asserter_type>, char>
//!\endcond
void read_n(receiver_type        && rcvr,
            input_rng_type       && input_rng,
            uint32_t       const    count,
            asserter_type        && asserter = detail::invocable_dummy)
{
    read_until(std::forward<receiver_type>(rcvr),
               std::forward<input_rng_type>(input_rng),
               [count = count] (auto) mutable { return (count-- > 0) ? false : true; },
               std::forward<asserter_type>(asserter));
}

// ----------------------------------------------------------------------------
// read_one
// ----------------------------------------------------------------------------

template <typename            receiver_type,
          input_range_concept input_rng_type,
          typename            asserter_type  = decltype(detail::invocable_dummy) &>
//!\cond
    requires (std::is_same_v<remove_cvref_t<receiver_type>, remove_cvref_t<decltype(std::ignore)>> ||
              output_iterator_concept<std::remove_reference_t<receiver_type>, char>) &&
              invocable_concept<std::remove_reference_t<asserter_type>, char>
//!\endcond
void read_one(receiver_type  && rcvr,
              input_rng_type && input_rng,
              asserter_type  && asserter = detail::invocable_dummy)
{
    read_n(std::forward<receiver_type>(rcvr),
           std::forward<input_rng_type>(input_rng),
           1,
           std::forward<asserter_type>(asserter));
}

//TODO(rrahn) add std::interface for ostream
// put
// write
} // namespace seqan3

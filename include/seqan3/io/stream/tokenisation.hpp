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

#include <seqan3/std/concept/core_language.hpp>
#include <seqan3/std/concept/iterator.hpp>
#include <seqan3/core/metafunction/basic.hpp>
#include <seqan3/core/metafunction/range.hpp>
#include <seqan3/io/exception.hpp>
#include <seqan3/range/concept.hpp>

namespace seqan3
{

// ----------------------------------------------------------------------------
// copy
// ----------------------------------------------------------------------------

/*!\brief      Copies elements from the input range to the output iterator, optionally checking certain conditions.
 * \ingroup    stream
 * \tparam     input_rng_type      The type of the input range; must satisfy seqan3::input_range_concept.
 * \tparam     out_iterator_type   The type of the output iterator; must satisfy seqan3::output_iterator_concept.
 * \tparam     stop_type           Type of the "stop" condition; must satisfy seqan3::predicate_concept (unless set to
 *                                 std::ignore).
 * \tparam     fail_type           Type of the "fail" condition; must satisfy seqan3::predicate_concept (unless set to
 *                                 std::ignore).
 * \tparam     skip_type           Type of the "skip" condition; must satisfy seqan3::predicate_concept (unless set to
 *                                 std::ignore).
 * \param[in]  input_range         The range to be parsed.
 * \param[out] output_it           The iterator to be written to.
 * \param[in]  stop_if             The condition on which to stop [optional]; usually a seqan3::parse_condition.
 * \param[in]  fail_if             The condition on which to fail [optional]; usually a seqan3::parse_condition.
 * \param[in]  skip_if             The condition on which to skip [optional]; usually a seqan3::parse_condition.
 * \throws     seqan3::parse_error If the fail condition is met or if the end of the input is reached without the
 *                                 stop condition being met.
 *
 * \details
 *
 * This function behaves similarly to std::copy, but it can perform certain condition checks to achieve fine-grained
 * tokenisation. The conditions are checked in the order: "stop", "fail", "skip". Any or all can be set to
 * std::ignore in which case the check will not be performed (this is the default if the parameters are not specified).
 *
 * ### Stop condition
 *
 * By default all elements of the input range are assigned to the output iterator, but if a "stop" condition is set,
 * it will instead stop if this condition is met. If the end of the input range is reached without satisfying
 * the condition an exception is thrown. To achieve a "stop condition is met or at end"-behaviour, you can make
 * the condition contain seqan3::is_end.
 *
 * ### Fail condition
 *
 * If you expect all elements read to satisfy certain conditions (e.g. be letters) and you want to raise an
 * exception if this is not the case, you can specify a fail condition. If it returns `true` on an element,
 * a seqan3::parse_error will be thrown.
 *
 * ### Skip condition
 *
 * If you wish to filter out certain elements without raising an exception, you can specify a skip condition.
 */
template <input_range_concept           input_rng_type,
          output_iterator_concept<char> out_iterator_type,
          typename                      stop_type = detail::ignore_t &,
          typename                      fail_type = detail::ignore_t &,
          typename                      skip_type = detail::ignore_t &>
//!\cond
    requires (predicate_concept<std::remove_reference_t<stop_type>, char> || decays_to_ignore_v<stop_type>) &&
             (predicate_concept<std::remove_reference_t<fail_type>, char> || decays_to_ignore_v<fail_type>) &&
             (predicate_concept<std::remove_reference_t<skip_type>, char> || decays_to_ignore_v<skip_type>)
//!\endcond
void copy(input_rng_type    && input_range,
          out_iterator_type && output_it,
          stop_type         && stop_if = std::ignore,
          fail_type         && fail_if = std::ignore,
          skip_type         && skip_if = std::ignore)
{
    //TODO extract is_eof from stop_type  to avoid redundant checks in the loop
    for (auto && input_char : input_range)
    {
        if constexpr (!decays_to_ignore_v<stop_type>)
            if (stop_if(input_char))
                return;

        if constexpr (!decays_to_ignore_v<fail_type>)
        {
            if (fail_if(input_char))
            {
                if constexpr (detail::parse_condition_concept<fail_type>)
                    throw parse_error{fail_type::msg.string()};
                else
                    throw parse_error{"Fail condition met while parsing character "s +
                                      detail::make_printable(input_char) + "."};
            }
        }

        if constexpr (!decays_to_ignore_v<skip_type>)
            if (skip_if(input_char))
                continue;

        *output_it = input_char;
        ++output_it;
    }

    //TODO and stop_type  does not contain is_eof
    if constexpr (!decays_to_ignore_v<stop_type>)
        throw unexpected_end_of_error{"Reached end of input while expecting more data."};
}

// ----------------------------------------------------------------------------
// copy_line
// ----------------------------------------------------------------------------

template <input_range_concept           input_rng_type,
          output_iterator_concept<char> out_iterator_type,
          typename                      fail_type = detail::ignore_t &,
          typename                      skip_type = detail::ignore_t &>
//!\cond
    requires (predicate_concept<std::remove_reference_t<fail_type>, char> || decays_to_ignore_v<fail_type>) &&
             (predicate_concept<std::remove_reference_t<skip_type>, char> || decays_to_ignore_v<skip_type>)
//!\endcond
void copy_line(input_rng_type    && input_range,
               out_iterator_type && output_it,
               fail_type         && fail_if = std::ignore,
               skip_type         && skip_if = std::ignore)
{
    copy(std::forward<input_rng_type>(input_range),
         std::forward<out_iterator_type>(output_it),
         is_char<'\n'>{} || is_char<'\r'>{},
         std::forward<fail_type>(fail_if),
         std::forward<skip_type>(skip_if));

    // Check if the statement was carriage return plus new-line.
    auto it = ranges::begin(input_rng);
    if (*it == '\r')
        if (*(++it) != '\n')  // consume the '\r' symbol and check if '\n' follows.
            throw parse_error{"Missing newline '\\n' character after reading '\\r' character."};
    ++it; // extract the newline character.
}

} // namespace seqan3

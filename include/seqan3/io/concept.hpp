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
 * \ingroup IO
 * \brief Stream concepts.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <type_traits>
#include <ios>

#include <range/v3/utility/associated_types.hpp>

namespace seqan3
{

/*!\addtogroup io
 * \{
 */

/*!\interface seqan3::ostream_concept <>
 * \brief Concept for output streams.
 *
 * An object is an output stream if it inherits from the [std::ios_base](http://en.cppreference.com/w/cpp/io/ios_base)
 * and supports the (un)formatted output function (\c operator<<) for a l-value of a given \c value_type.
 * It further needs to define the public member types as described in the [STD](http://en.cppreference.com/w/cpp/io/basic_ios).
 */
//!\cond
template <typename stream_type, typename value_type>
concept bool ostream_concept = std::is_base_of_v<std::ios_base, stream_type> &&
                               requires (stream_type & os, value_type & val)
{
    typename stream_type::char_type;
    typename stream_type::traits_type;
    typename stream_type::int_type;
    typename stream_type::pos_type;
    typename stream_type::off_type;

    { os << val } -> std::basic_ostream<typename stream_type::char_type, typename stream_type::traits_type> &;
};
//!\endcond

/*!\interface seqan3::istream_concept <>
 * \brief Concept for input streams.
 *
 * An object is an input stream if it inherits from the [std::ios_base](http://en.cppreference.com/w/cpp/io/ios_base)
 * and supports the (un)formatted input function (\c operator>>) for a l-value of a given \c value_type.
 * It further needs to define the public member types as described in the [STD](http://en.cppreference.com/w/cpp/io/basic_ios).
 */
//!\cond
template <typename stream_type, typename value_type>
concept bool istream_concept = std::is_base_of_v<std::ios_base, stream_type> &&
                               requires (stream_type & os, value_type & val)
{
    typename stream_type::char_type;
    typename stream_type::traits_type;
    typename stream_type::int_type;
    typename stream_type::pos_type;
    typename stream_type::off_type;

    { os >> val } -> std::basic_istream<typename stream_type::char_type, typename stream_type::traits_type> &;
};
//!\endcond

/*!\interface seqan3::stream_concept <>
 * \extends seqan3::istream_concept
 * \extends seqan3::ostream_concept
 * \brief Concept for i/o streams permitting both directions .
 *
 * An object fulfilling the requirements of a stream concept can be used to stream in both (input and output)
 * directions.
 */
//!\cond
template <typename stream_type, typename value_type>
concept bool stream_concept = ostream_concept<stream_type, value_type> &&
                              istream_concept<stream_type, value_type>;
//!\endcond

//!\}

} // namespace seqan3

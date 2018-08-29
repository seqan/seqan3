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
 * \brief Stream concepts.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <ios>
#include <type_traits>

#include <range/v3/utility/associated_types.hpp>

#include <seqan3/core/platform.hpp>

namespace seqan3
{

/*!\addtogroup io
 * \{
 */

/*!\interface seqan3::ostream_concept <>
 * \brief Concept for output streams.
 *
 * An object is an output stream if it inherits from the [std::ios_base](http://en.cppreference.com/w/cpp/io/ios_base)
 * and supports the (un)formatted output function (`operator<<``) for a l-value of a given `value_type`.
 * It further needs to define the public member types as described in the [STD](http://en.cppreference.com/w/cpp/io/basic_ios).
 */
//!\cond
template <typename stream_type, typename value_type>
concept ostream_concept = std::is_base_of_v<std::ios_base, stream_type> &&
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

/*!\name Requirements for seqan3::ostream_concept
 * \relates seqan3::ostream_concept
 * \brief You can expect these member types and the free function on all types that satisfy seqan3::ostream_concept.
 * \{
 */
/*!\fn      std::basic_ostream<char_type, traits_type> & operator<<(value_type val);
 * \brief   (un)-formatted output operator for the respective type on the underlying stream.
 * \ingroup io
 * \param   val The value to write into the stream.
 * \returns A reference to a std::basic_ostream<char_type, traits_type>.
 *
 * \details
 *
 * The `char_type` and `traits_type` are inferred from the given `ostream`.
 *
 * \attention This is a concept requirement, not an actual function (however types satisfying this concept
 * will provide an implementation).
 */

 /*!\typedef typename stream::char_type char_type
  * \memberof seqan3::ostream_concept
  * \brief Declares the associated char type.
  */

 /*!\typedef typename stream::traits_type traits_type
  * \memberof seqan3::ostream_concept
  * \brief Declares the associated traits type.
  */

 /*!\typedef typename stream::int_type int_type
  * \memberof seqan3::ostream_concept
  * \brief Declares the associated int type.
  */

 /*!\typedef typename stream::pos_type pos_type
  * \memberof seqan3::ostream_concept
  * \brief Declares the associated pos type.
  */

 /*!\typedef typename stream::off_type off_type
  * \memberof seqan3::ostream_concept
  * \brief Declares the associated off type.
  */
//!\}

/*!\interface seqan3::istream_concept <>
 * \brief Concept for input streams.
 *
 * An object is an input stream if it inherits from the [std::ios_base](http://en.cppreference.com/w/cpp/io/ios_base)
 * and supports the (un)formatted input function (`operator>>`) for a l-value of a given `value_type`.
 * It further needs to define the public member types as described in the [STD](http://en.cppreference.com/w/cpp/io/basic_ios).
 */
//!\cond
template <typename stream_type, typename value_type>
concept istream_concept = std::is_base_of_v<std::ios_base, stream_type> &&
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

/*!\name Requirements for seqan3::istream_concept
 * \relates seqan3::istream_concept
 * \brief You can expect these member types and the free function on all types that satisfy seqan3::istream_concept.
 * \{
 */
/*!\fn      std::basic_istream<char_type, traits_type> & operator>>(value_type val);
 * \brief   (un)-formatted input operator for the respective type on the underlying stream.
 * \ingroup io
 * \param   val The value to read from the stream.
 * \returns A reference to a std::basic_istream<char_type, traits_type>.
 *
 * \details
 *
 * The `char_type` and `traits_type` are inferred from the given `istream`.
 *
 * \attention This is a concept requirement, not an actual function (however types satisfying this concept
 * will provide an implementation).
 */

 /*!\typedef typename stream::char_type char_type
  * \memberof seqan3::istream_concept
  * \brief Declares the associated char type.
  */

 /*!\typedef typename stream::traits_type traits_type
  * \memberof seqan3::istream_concept
  * \brief Declares the associated traits type.
  */

 /*!\typedef typename stream::int_type int_type
  * \memberof seqan3::istream_concept
  * \brief Declares the associated int type.
  */

 /*!\typedef typename stream::pos_type pos_type
  * \memberof seqan3::istream_concept
  * \brief Declares the associated pos type.
  */

 /*!\typedef typename stream::off_type off_type
  * \memberof seqan3::istream_concept
  * \brief Declares the associated off type.
  */
//!\}

/*!\interface seqan3::stream_concept <>
 * \extends seqan3::istream_concept
 * \extends seqan3::ostream_concept
 * \brief Concept for i/o streams permitting both directions.
 *
 * An object satisfying the requirements of a stream concept can be used to stream in both (input and output)
 * directions.
 */
//!\cond
template <typename stream_type, typename value_type>
concept stream_concept = ostream_concept<stream_type, value_type> &&
                              istream_concept<stream_type, value_type>;
//!\endcond

//!\}

} // namespace seqan3

// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

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
/*!\interface seqan3::ostream_concept <>
 * \ingroup stream
 * \brief Concept for output streams.
 *
 * An object is an output stream if it inherits from the [std::ios_base](http://en.cppreference.com/w/cpp/io/ios_base)
 * and supports the (un)formatted output function (`operator<<``) for a l-value of a given `value_type`.
 * It further needs to define the public member types as described in the [STD](http://en.cppreference.com/w/cpp/io/basic_ios).
 */
//!\cond
template <typename stream_type, typename value_type>
SEQAN3_CONCEPT ostream_concept = std::is_base_of_v<std::ios_base, std::remove_reference_t<stream_type>> &&
                               requires (stream_type & os, value_type & val)
{
    typename std::remove_reference_t<stream_type>::char_type;
    typename std::remove_reference_t<stream_type>::traits_type;
    typename std::remove_reference_t<stream_type>::int_type;
    typename std::remove_reference_t<stream_type>::pos_type;
    typename std::remove_reference_t<stream_type>::off_type;

    { os << val } -> std::basic_ostream<typename std::remove_reference_t<stream_type>::char_type,
                                        typename std::remove_reference_t<stream_type>::traits_type> &;
};

template <typename stream_type>
SEQAN3_CONCEPT ostream_concept2 = requires { typename std::remove_reference_t<stream_type>::char_type; } &&
                           ostream_concept<stream_type, typename std::remove_reference_t<stream_type>::char_type>;
//!\endcond

/*!\name Requirements for seqan3::ostream_concept
 * \relates seqan3::ostream_concept
 * \brief You can expect these member types and the free function on all types that satisfy seqan3::ostream_concept.
 * \{
 */
/*!\fn      std::basic_ostream<char_type, traits_type> & operator<<(value_type val);
 * \brief   (un)-formatted output operator for the respective type on the underlying stream.
 * \relates seqan3::ostream_concept
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
 * \ingroup stream
 * \brief Concept for input streams.
 *
 * An object is an input stream if it inherits from the [std::ios_base](http://en.cppreference.com/w/cpp/io/ios_base)
 * and supports the (un)formatted input function (`operator>>`) for a l-value of a given `value_type`.
 * It further needs to define the public member types as described in the [STD](http://en.cppreference.com/w/cpp/io/basic_ios).
 */
//!\cond
template <typename stream_type, typename value_type>
SEQAN3_CONCEPT istream_concept = std::is_base_of_v<std::ios_base, std::remove_reference_t<stream_type>> &&
                               requires (stream_type & is, value_type & val)
{
    typename std::remove_reference_t<stream_type>::char_type;
    typename std::remove_reference_t<stream_type>::traits_type;
    typename std::remove_reference_t<stream_type>::int_type;
    typename std::remove_reference_t<stream_type>::pos_type;
    typename std::remove_reference_t<stream_type>::off_type;

    { is >> val } -> std::basic_istream<typename std::remove_reference_t<stream_type>::char_type,
                                        typename std::remove_reference_t<stream_type>::traits_type> &;
};

template <typename stream_type>
SEQAN3_CONCEPT istream_concept2 = requires { typename std::remove_reference_t<stream_type>::char_type; } &&
                           istream_concept<stream_type, typename std::remove_reference_t<stream_type>::char_type>;
//!\endcond

/*!\name Requirements for seqan3::istream_concept
 * \relates seqan3::istream_concept
 * \brief You can expect these member types and the free function on all types that satisfy seqan3::istream_concept.
 * \{
 */
/*!\fn      std::basic_istream<char_type, traits_type> & operator>>(value_type val);
 * \brief   (un)-formatted input operator for the respective type on the underlying stream.
 * \relates seqan3::istream_concept
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
 * \ingroup stream
 *
 * An object satisfying the requirements of a stream concept can be used to stream in both (input and output)
 * directions.
 */
//!\cond
template <typename stream_type, typename value_type>
SEQAN3_CONCEPT stream_concept = ostream_concept<stream_type, value_type> &&
                              istream_concept<stream_type, value_type>;

//!\endcond

} // namespace seqan3

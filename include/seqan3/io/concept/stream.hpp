// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
//
// Copyright (c) 2006-2017, Knut Reinert, FU Berlin
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
// ==========================================================================

#pragma once

/*! \file io/concept/strem_concept.hpp
 *  \ingroup io
 *  \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

 #include <fstream>
 #include <sstream>

 namespace seqan3::detail
 {
 // basic meta_functions:
 template <typename t>
 struct stream_traits
 {
     using char_type   = void;
     using traits_type = void;
     using int_type    = void;
     using off_type    = void;
     using pos_type    = void;
 };

 template <typename char_t, typename traits_t>
 struct stream_traits<std::basic_iostream<char_t, traits_t>>
 {
     using stream_ = std::basic_iostream<char_t, traits_t>;
     using char_type   = typename stream_::char_type;
     using traits_type = typename stream_::traits_type;
     using int_type    = typename stream_::int_type;
     using off_type    = typename stream_::off_type;
     using pos_type    = typename stream_::pos_type;
 };

 template <typename char_t, typename traits_t>
 struct stream_traits<std::basic_fstream<char_t, traits_t>>
 {
     using stream_ = std::basic_fstream<char_t, traits_t>;

     using char_type   = typename stream_::char_type;
     using traits_type = typename stream_::traits_type;
     using int_type    = typename stream_::int_type;
     using off_type    = typename stream_::off_type;
     using pos_type    = typename stream_::pos_type;
 };

 template <typename char_t, typename traits_t>
 struct stream_traits<std::basic_stringstream<char_t, traits_t>>
 {
     using stream_ = std::basic_stringstream<char_t, traits_t>;

     using char_type   = typename stream_::char_type;
     using traits_type = typename stream_::traits_type;
     using int_type    = typename stream_::int_type;
     using off_type    = typename stream_::off_type;
     using pos_type    = typename stream_::pos_type;
 };

 template <typename char_t, typename traits_t>
 struct stream_traits<std::basic_istream<char_t, traits_t>>
 {
     using stream_ = std::basic_istream<char_t, traits_t>;

     using char_type   = typename stream_::char_type;
     using traits_type = typename stream_::traits_type;
     using int_type    = typename stream_::int_type;
     using off_type    = typename stream_::off_type;
     using pos_type    = typename stream_::pos_type;
 };

 template <typename char_t, typename traits_t>
 struct stream_traits<std::basic_ifstream<char_t, traits_t>>
 {
     using stream_ = std::basic_ifstream<char_t, traits_t>;

     using char_type   = typename stream_::char_type;
     using traits_type = typename stream_::traits_type;
     using int_type    = typename stream_::int_type;
     using off_type    = typename stream_::off_type;
     using pos_type    = typename stream_::pos_type;
 };

 template <typename char_t, typename traits_t>
 struct stream_traits<std::basic_istringstream<char_t, traits_t>>
 {
     using stream_ = std::basic_istringstream<char_t, traits_t>;

     using char_type   = typename stream_::char_type;
     using traits_type = typename stream_::traits_type;
     using int_type    = typename stream_::int_type;
     using off_type    = typename stream_::off_type;
     using pos_type    = typename stream_::pos_type;
 };

 template <typename char_t, typename traits_t>
 struct stream_traits<std::basic_ostream<char_t, traits_t>>
 {
     using stream_ = std::basic_ostream<char_t, traits_t>;

     using char_type   = typename stream_::char_type;
     using traits_type = typename stream_::traits_type;
     using int_type    = typename stream_::int_type;
     using off_type    = typename stream_::off_type;
     using pos_type    = typename stream_::pos_type;
 };

 template <typename char_t, typename traits_t>
 struct stream_traits<std::basic_ofstream<char_t, traits_t>>
 {
     using stream_ = std::basic_ofstream<char_t, traits_t>;

     using char_type   = typename stream_::char_type;
     using traits_type = typename stream_::traits_type;
     using int_type    = typename stream_::int_type;
     using off_type    = typename stream_::off_type;
     using pos_type    = typename stream_::pos_type;
 };

 template <typename char_t, typename traits_t>
 struct stream_traits<std::basic_ostringstream<char_t, traits_t>>
 {
     using stream_ = std::basic_ostringstream<char_t, traits_t>;
     using char_type   = typename stream_::char_type;
     using traits_type = typename stream_::traits_type;
     using int_type    = typename stream_::int_type;
     using off_type    = typename stream_::off_type;
     using pos_type    = typename stream_::pos_type;
 };
 }  // namespace seqan3::detail

namespace seqan3
{

template <typename stream_t>
struct char_type
{
    using type = typename detail::stream_traits<std::decay_t<stream_t>>::char_type;
};

template <typename stream_t>
using char_type_t =  typename char_type<stream_t>::type;

template <typename stream_t>
struct traits_type
{
    using type = typename detail::stream_traits<std::decay_t<stream_t>>::traits_type;
};

template <typename stream_t>
using traits_type_t =  typename traits_type<stream_t>::type;

template <typename stream_t>
struct int_type
{
    using type = typename detail::stream_traits<std::decay_t<stream_t>>::int_type;
};

template <typename stream_t>
using int_type_t =  typename int_type<stream_t>::type;

template <typename stream_t>
struct off_type
{
    using type = typename detail::stream_traits<std::decay_t<stream_t>>::off_type;
};

template <typename stream_t>
using off_type_t =  typename off_type<stream_t>::type;

template <typename stream_t>
struct pos_type
{
    using type = typename detail::stream_traits<std::decay_t<stream_t>>::pos_type;
};

template <typename stream_t>
using pos_type_t =  typename pos_type<stream_t>::type;

template <typename t>
concept bool input_stream_concept =
    requires(t stream, typename t::char_type val)
{
    typename char_type<t>::type;
    typename traits_type<t>::type;
    typename int_type<t>::type;
    typename pos_type<t>::type;
    typename off_type<t>::type;

    // formatted input.
    { operator>>(stream, val) } -> auto&;
};

template <typename t>
concept bool output_stream_concept =
    requires(t stream)
{
    typename char_type<t>::type;
    typename traits_type<t>::type;
    typename int_type<t>::type;
    typename pos_type<t>::type;
    typename off_type<t>::type;

    // formatted output.
    { operator<<(stream, typename t::char_type{}) } -> auto&;
};

template <typename t>
concept bool bidirectional_stream_concept = input_stream_concept<t> &&
                                            output_stream_concept<t>;

}  // namespace seqan3

//#ifndef NDEBUG

namespace seqan3::detail
{
#include <type_traits>
#include <sstream>
#include <vector>

// check input_stream_concept
static_assert(input_stream_concept<std::istringstream>);
static_assert(!input_stream_concept<std::ostringstream>);
static_assert(input_stream_concept<std::stringstream>);
static_assert(!input_stream_concept<std::vector<char>>);

// check output_stream_concept
static_assert(!output_stream_concept<std::istringstream>);
static_assert(output_stream_concept<std::ostringstream>);
static_assert(output_stream_concept<std::stringstream>);
static_assert(!output_stream_concept<std::vector<char>>);

// check bidirectional_stream_concept
static_assert(!bidirectional_stream_concept<std::istringstream>);
static_assert(!bidirectional_stream_concept<std::ostringstream>);
static_assert(bidirectional_stream_concept<std::stringstream>);
static_assert(!bidirectional_stream_concept<std::vector<char>>);
}

//#endif // NDEBUG

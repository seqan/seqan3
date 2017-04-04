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
// Author: Rene Rahn <rene.rahn AT fu-berlin.de>
// ============================================================================

#pragma once

#ifndef NDEBUG

#include <iterator>
#include <list>
#include <forward_list>
#include <vector>

#include <seqan3/core/concepts/iterator.hpp>

/* Checks the iterator concepts */

namespace seqan3::detail
{
using check_concept_for_input_iterator               = std::istream_iterator<char>;
using check_concept_for_output_iterator              = ranges::ostream_iterator<char>;
using check_concept_for_forward_iterator             = std::forward_list<char>::iterator;
using check_concept_for_bidirectional_iterator       = std::list<char>::iterator;
using check_concept_for_random_access_iterator       = std::vector<char>::iterator;

using check_concept_for_forward_iterator_const       = std::forward_list<char>::const_iterator;
using check_concept_for_bidirectional_iterator_const = std::list<char>::const_iterator;
using check_concept_for_random_access_iterator_const = std::vector<char>::const_iterator;

// Weakly weakly_incrementable, semi-regular, weakly_equality_comparable
template <typename value_t>
struct test_sentinel
{
    using value_type      = value_t;
    using difference_type = size_t;

    value_type val{};
};

template <typename iterator_t, typename value_t>
inline bool operator==(iterator_t const & i,
                       test_sentinel<value_t> const & s)
{
    return *i == s.val;
}

template <typename iterator_t, typename value_t>
inline bool operator!=(iterator_t const & i,
                       test_sentinel<value_t> const & s)
{
    return !(i == s);
}

template <typename value_t, typename iterator_t>
inline bool operator==(test_sentinel<value_t> const & s,
                       iterator_t const & i)
{
    return *i == s.val;
}

template <typename value_t, typename iterator_t>
inline bool operator!=(test_sentinel<value_t> const & s,
                       iterator_t const & i)
{
    return !(i == s);
}

template <typename iterator_type>
struct value
{
    using type = typename std::iterator_traits<iterator_type>::value_type;
};

template <typename value_t, typename ...ts>
struct value<ranges::ostream_iterator<value_t, ts...>>
{
    using type = value_t;
};

template <typename iterator_type>
using value_type_t = typename value<iterator_type>::type;

template <typename iterator_type>
struct test_sized_sentinel : public test_sentinel<value_type_t<iterator_type>>
{
    using difference_type = typename iterator_type::difference_type;

    iterator_type _pos;
};

template <typename iterator_t>
    requires random_access_iterator_concept<iterator_t>
inline typename test_sized_sentinel<iterator_t>::difference_type
operator-(test_sized_sentinel<iterator_t> const & s,
          iterator_t const & i)
{
    return s._pos - i;
}

template <typename iterator_t>
    requires random_access_iterator_concept<iterator_t>
inline typename test_sized_sentinel<iterator_t>::difference_type
operator-(iterator_t const & i,
          test_sized_sentinel<iterator_t> const & s)
{
    return i - s._pos;
}
}  // namspace seqan3::detail

// Check readable_concept
static_assert(seqan3::readable_concept<seqan3::detail::check_concept_for_input_iterator>);
static_assert(!seqan3::readable_concept<seqan3::detail::check_concept_for_output_iterator>);
static_assert(seqan3::readable_concept<seqan3::detail::check_concept_for_forward_iterator>);
static_assert(seqan3::readable_concept<seqan3::detail::check_concept_for_bidirectional_iterator>);
static_assert(seqan3::readable_concept<seqan3::detail::check_concept_for_random_access_iterator>);
static_assert(seqan3::readable_concept<seqan3::detail::check_concept_for_forward_iterator_const>);
static_assert(seqan3::readable_concept<seqan3::detail::check_concept_for_bidirectional_iterator_const>);
static_assert(seqan3::readable_concept<seqan3::detail::check_concept_for_random_access_iterator_const>);

// Check writable_concept
static_assert(!seqan3::writable_concept<seqan3::detail::check_concept_for_input_iterator, char>);
static_assert(seqan3::writable_concept<seqan3::detail::check_concept_for_output_iterator, char>);
static_assert(seqan3::writable_concept<seqan3::detail::check_concept_for_forward_iterator, char>);
static_assert(!seqan3::writable_concept<seqan3::detail::check_concept_for_forward_iterator_const, char>);
static_assert(seqan3::writable_concept<seqan3::detail::check_concept_for_bidirectional_iterator, char>);
static_assert(!seqan3::writable_concept<seqan3::detail::check_concept_for_bidirectional_iterator_const, char>);
static_assert(seqan3::writable_concept<seqan3::detail::check_concept_for_random_access_iterator, char>);
static_assert(!seqan3::writable_concept<seqan3::detail::check_concept_for_random_access_iterator_const, char>);

// Check weakly_incrementable_concept
static_assert(seqan3::weakly_incrementable_concept<seqan3::detail::check_concept_for_input_iterator>);
static_assert(seqan3::weakly_incrementable_concept<seqan3::detail::check_concept_for_output_iterator>);
static_assert(seqan3::weakly_incrementable_concept<seqan3::detail::check_concept_for_forward_iterator>);
static_assert(seqan3::weakly_incrementable_concept<seqan3::detail::check_concept_for_forward_iterator_const>);
static_assert(seqan3::weakly_incrementable_concept<seqan3::detail::check_concept_for_bidirectional_iterator>);
static_assert(seqan3::weakly_incrementable_concept<seqan3::detail::check_concept_for_bidirectional_iterator_const>);
static_assert(seqan3::weakly_incrementable_concept<seqan3::detail::check_concept_for_random_access_iterator>);
static_assert(seqan3::weakly_incrementable_concept<seqan3::detail::check_concept_for_random_access_iterator_const>);

// Check incrementable_concept
static_assert(seqan3::incrementable_concept<seqan3::detail::check_concept_for_input_iterator>);
static_assert(!seqan3::incrementable_concept<seqan3::detail::check_concept_for_output_iterator>);
static_assert(seqan3::incrementable_concept<seqan3::detail::check_concept_for_forward_iterator>);
static_assert(seqan3::incrementable_concept<seqan3::detail::check_concept_for_forward_iterator_const>);
static_assert(seqan3::incrementable_concept<seqan3::detail::check_concept_for_bidirectional_iterator>);
static_assert(seqan3::incrementable_concept<seqan3::detail::check_concept_for_bidirectional_iterator_const>);
static_assert(seqan3::incrementable_concept<seqan3::detail::check_concept_for_random_access_iterator>);
static_assert(seqan3::incrementable_concept<seqan3::detail::check_concept_for_random_access_iterator_const>);

// Check iterator_concept
static_assert(seqan3::iterator_concept<seqan3::detail::check_concept_for_input_iterator>);
static_assert(seqan3::iterator_concept<seqan3::detail::check_concept_for_output_iterator>);
static_assert(seqan3::iterator_concept<seqan3::detail::check_concept_for_forward_iterator>);
static_assert(seqan3::iterator_concept<seqan3::detail::check_concept_for_forward_iterator_const>);
static_assert(seqan3::iterator_concept<seqan3::detail::check_concept_for_bidirectional_iterator>);
static_assert(seqan3::iterator_concept<seqan3::detail::check_concept_for_bidirectional_iterator_const>);
static_assert(seqan3::iterator_concept<seqan3::detail::check_concept_for_random_access_iterator>);
static_assert(seqan3::iterator_concept<seqan3::detail::check_concept_for_random_access_iterator_const>);

// Check sentinel_concept
static_assert(seqan3::sentinel_concept<seqan3::detail::test_sentinel<char>,
                                       seqan3::detail::check_concept_for_input_iterator>);
static_assert(seqan3::sentinel_concept<seqan3::detail::test_sentinel<char>,
                                       seqan3::detail::check_concept_for_output_iterator>);
static_assert(seqan3::sentinel_concept<seqan3::detail::test_sentinel<char>,
                                       seqan3::detail::check_concept_for_forward_iterator>);
static_assert(seqan3::sentinel_concept<seqan3::detail::test_sentinel<char>,
                                       seqan3::detail::check_concept_for_forward_iterator_const>);
static_assert(seqan3::sentinel_concept<seqan3::detail::test_sentinel<char>,
                                       seqan3::detail::check_concept_for_bidirectional_iterator>);
static_assert(seqan3::sentinel_concept<seqan3::detail::test_sentinel<char>,
                                       seqan3::detail::check_concept_for_bidirectional_iterator_const>);
static_assert(seqan3::sentinel_concept<seqan3::detail::test_sentinel<char>,
                                       seqan3::detail::check_concept_for_random_access_iterator>);
static_assert(seqan3::sentinel_concept<seqan3::detail::test_sentinel<char>,
                                       seqan3::detail::check_concept_for_random_access_iterator_const>);

static_assert(seqan3::sentinel_concept<seqan3::detail::test_sized_sentinel<seqan3::detail::check_concept_for_input_iterator>,
                                       seqan3::detail::check_concept_for_input_iterator>);
static_assert(std::is_same_v<char, seqan3::detail::value_type_t<ranges::ostream_iterator<char>>>);
static_assert(seqan3::sentinel_concept<seqan3::detail::test_sized_sentinel<seqan3::detail::check_concept_for_output_iterator>,
                                       seqan3::detail::check_concept_for_output_iterator>);
static_assert(seqan3::sentinel_concept<seqan3::detail::test_sized_sentinel<seqan3::detail::check_concept_for_forward_iterator>,
                                       seqan3::detail::check_concept_for_forward_iterator>);
static_assert(seqan3::sentinel_concept<seqan3::detail::test_sized_sentinel<seqan3::detail::check_concept_for_forward_iterator_const>,
                                       seqan3::detail::check_concept_for_forward_iterator_const>);
static_assert(seqan3::sentinel_concept<seqan3::detail::test_sized_sentinel<seqan3::detail::check_concept_for_bidirectional_iterator>,
                                       seqan3::detail::check_concept_for_bidirectional_iterator>);
static_assert(seqan3::sentinel_concept<seqan3::detail::test_sized_sentinel<seqan3::detail::check_concept_for_bidirectional_iterator_const>,
                                       seqan3::detail::check_concept_for_bidirectional_iterator_const>);
static_assert(seqan3::sentinel_concept<seqan3::detail::test_sized_sentinel<seqan3::detail::check_concept_for_random_access_iterator>,
                                       seqan3::detail::check_concept_for_random_access_iterator>);
static_assert(seqan3::sentinel_concept<seqan3::detail::test_sized_sentinel<seqan3::detail::check_concept_for_random_access_iterator_const>,
                                       seqan3::detail::check_concept_for_random_access_iterator_const>);

// Check sized_sentinel_concept
static_assert(!seqan3::sized_sentinel_concept<seqan3::detail::test_sentinel<char>,
                                              seqan3::detail::check_concept_for_input_iterator>);
static_assert(!seqan3::sized_sentinel_concept<seqan3::detail::test_sentinel<char>,
                                              seqan3::detail::check_concept_for_output_iterator>);
static_assert(!seqan3::sized_sentinel_concept<seqan3::detail::test_sentinel<char>,
                                              seqan3::detail::check_concept_for_forward_iterator>);
static_assert(!seqan3::sized_sentinel_concept<seqan3::detail::test_sentinel<char>,
                                              seqan3::detail::check_concept_for_forward_iterator_const>);
static_assert(!seqan3::sized_sentinel_concept<seqan3::detail::test_sentinel<char>,
                                              seqan3::detail::check_concept_for_bidirectional_iterator>);
static_assert(!seqan3::sized_sentinel_concept<seqan3::detail::test_sentinel<char>,
                                              seqan3::detail::check_concept_for_bidirectional_iterator_const>);
static_assert(!seqan3::sized_sentinel_concept<seqan3::detail::test_sentinel<char>,
                                              seqan3::detail::check_concept_for_random_access_iterator>);
static_assert(!seqan3::sized_sentinel_concept<seqan3::detail::test_sentinel<char>,
                                              seqan3::detail::check_concept_for_random_access_iterator_const>);

static_assert(!seqan3::sized_sentinel_concept<seqan3::detail::test_sized_sentinel<seqan3::detail::check_concept_for_input_iterator>,
                                             seqan3::detail::check_concept_for_input_iterator>);
static_assert(!seqan3::sized_sentinel_concept<seqan3::detail::test_sized_sentinel<seqan3::detail::check_concept_for_output_iterator>,
                                             seqan3::detail::check_concept_for_output_iterator>);
static_assert(!seqan3::sized_sentinel_concept<seqan3::detail::test_sized_sentinel<seqan3::detail::check_concept_for_forward_iterator>,
                                             seqan3::detail::check_concept_for_forward_iterator>);
static_assert(!seqan3::sized_sentinel_concept<seqan3::detail::test_sized_sentinel<seqan3::detail::check_concept_for_forward_iterator_const>,
                                             seqan3::detail::check_concept_for_forward_iterator_const>);
static_assert(!seqan3::sized_sentinel_concept<seqan3::detail::test_sized_sentinel<seqan3::detail::check_concept_for_bidirectional_iterator>,
                                             seqan3::detail::check_concept_for_bidirectional_iterator>);
static_assert(!seqan3::sized_sentinel_concept<seqan3::detail::test_sized_sentinel<seqan3::detail::check_concept_for_bidirectional_iterator_const>,
                                             seqan3::detail::check_concept_for_bidirectional_iterator_const>);
static_assert(seqan3::sized_sentinel_concept<seqan3::detail::test_sized_sentinel<seqan3::detail::check_concept_for_random_access_iterator>,
                                             seqan3::detail::check_concept_for_random_access_iterator>);
static_assert(seqan3::sized_sentinel_concept<seqan3::detail::test_sized_sentinel<seqan3::detail::check_concept_for_random_access_iterator_const>,
                                             seqan3::detail::check_concept_for_random_access_iterator_const>);

// Check output_iterator_concept
static_assert(!seqan3::output_iterator_concept<seqan3::detail::check_concept_for_input_iterator, char>);
static_assert(seqan3::output_iterator_concept<seqan3::detail::check_concept_for_output_iterator, char>);
static_assert(seqan3::output_iterator_concept<seqan3::detail::check_concept_for_forward_iterator, char>);
static_assert(!seqan3::output_iterator_concept<seqan3::detail::check_concept_for_forward_iterator_const, char>);
static_assert(seqan3::output_iterator_concept<seqan3::detail::check_concept_for_bidirectional_iterator, char>);
static_assert(!seqan3::output_iterator_concept<seqan3::detail::check_concept_for_bidirectional_iterator_const, char>);
static_assert(seqan3::output_iterator_concept<seqan3::detail::check_concept_for_random_access_iterator, char>);
static_assert(!seqan3::output_iterator_concept<seqan3::detail::check_concept_for_random_access_iterator_const, char>);

// Check input_iterator_concept
static_assert(seqan3::input_iterator_concept<seqan3::detail::check_concept_for_input_iterator>);
static_assert(!seqan3::input_iterator_concept<seqan3::detail::check_concept_for_output_iterator>);
static_assert(seqan3::input_iterator_concept<seqan3::detail::check_concept_for_forward_iterator>);
static_assert(seqan3::input_iterator_concept<seqan3::detail::check_concept_for_forward_iterator_const>);
static_assert(seqan3::input_iterator_concept<seqan3::detail::check_concept_for_bidirectional_iterator>);
static_assert(seqan3::input_iterator_concept<seqan3::detail::check_concept_for_bidirectional_iterator_const>);
static_assert(seqan3::input_iterator_concept<seqan3::detail::check_concept_for_random_access_iterator>);
static_assert(seqan3::input_iterator_concept<seqan3::detail::check_concept_for_random_access_iterator_const>);

// Check forward_iterator_concept
static_assert(!seqan3::forward_iterator_concept<seqan3::detail::check_concept_for_input_iterator>);
static_assert(!seqan3::forward_iterator_concept<seqan3::detail::check_concept_for_output_iterator>);
static_assert(seqan3::forward_iterator_concept<seqan3::detail::check_concept_for_forward_iterator>);
static_assert(seqan3::forward_iterator_concept<seqan3::detail::check_concept_for_forward_iterator_const>);
static_assert(seqan3::forward_iterator_concept<seqan3::detail::check_concept_for_bidirectional_iterator>);
static_assert(seqan3::forward_iterator_concept<seqan3::detail::check_concept_for_bidirectional_iterator_const>);
static_assert(seqan3::forward_iterator_concept<seqan3::detail::check_concept_for_random_access_iterator>);
static_assert(seqan3::forward_iterator_concept<seqan3::detail::check_concept_for_random_access_iterator_const>);

// Check bidirectional_iterator_concept
static_assert(!seqan3::bidirectional_iterator_concept<seqan3::detail::check_concept_for_input_iterator>);
static_assert(!seqan3::bidirectional_iterator_concept<seqan3::detail::check_concept_for_output_iterator>);
static_assert(!seqan3::bidirectional_iterator_concept<seqan3::detail::check_concept_for_forward_iterator>);
static_assert(!seqan3::bidirectional_iterator_concept<seqan3::detail::check_concept_for_forward_iterator_const>);
static_assert(seqan3::bidirectional_iterator_concept<seqan3::detail::check_concept_for_bidirectional_iterator>);
static_assert(seqan3::bidirectional_iterator_concept<seqan3::detail::check_concept_for_bidirectional_iterator_const>);
static_assert(seqan3::bidirectional_iterator_concept<seqan3::detail::check_concept_for_random_access_iterator>);
static_assert(seqan3::bidirectional_iterator_concept<seqan3::detail::check_concept_for_random_access_iterator_const>);

// Check random_access_iterator_concept
static_assert(!seqan3::random_access_iterator_concept<seqan3::detail::check_concept_for_input_iterator>);
static_assert(!seqan3::random_access_iterator_concept<seqan3::detail::check_concept_for_output_iterator>);
static_assert(!seqan3::random_access_iterator_concept<seqan3::detail::check_concept_for_forward_iterator>);
static_assert(!seqan3::random_access_iterator_concept<seqan3::detail::check_concept_for_forward_iterator_const>);
static_assert(!seqan3::random_access_iterator_concept<seqan3::detail::check_concept_for_bidirectional_iterator>);
static_assert(!seqan3::random_access_iterator_concept<seqan3::detail::check_concept_for_bidirectional_iterator_const>);
static_assert(seqan3::random_access_iterator_concept<seqan3::detail::check_concept_for_random_access_iterator>);
static_assert(seqan3::random_access_iterator_concept<seqan3::detail::check_concept_for_random_access_iterator_const>);

#endif // NDEBUG

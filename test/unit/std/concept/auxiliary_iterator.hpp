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

#pragma once

#include <iterator>
#include <list>
#include <forward_list>
#include <vector>

#include <seqan3/std/iterator>
#include <seqan3/std/ranges>

using input_iterator               = std::istream_iterator<char>;
using output_iterator              = std::ranges::ostream_iterator<char>;
using forward_iterator             = std::forward_list<char>::iterator;
using bidirectional_iterator       = std::list<char>::iterator;
using random_access_iterator       = std::vector<char>::iterator;
using forward_iterator_const       = std::forward_list<char>::const_iterator;
using bidirectional_iterator_const = std::list<char>::const_iterator;
using random_access_iterator_const = std::vector<char>::const_iterator;

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
struct value<std::ranges::ostream_iterator<value_t, ts...>>
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
    requires std::RandomAccessIterator<iterator_t>
inline typename test_sized_sentinel<iterator_t>::difference_type
operator-(test_sized_sentinel<iterator_t> const & s,
          iterator_t const & i)
{
    return s._pos - i;
}

template <typename iterator_t>
    requires std::RandomAccessIterator<iterator_t>
inline typename test_sized_sentinel<iterator_t>::difference_type
operator-(iterator_t const & i,
          test_sized_sentinel<iterator_t> const & s)
{
    return i - s._pos;
}

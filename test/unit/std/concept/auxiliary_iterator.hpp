// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <forward_list>
#include <seqan3/std/iterator>
#include <list>
#include <vector>

using input_iterator = std::istream_iterator<char>;
using output_iterator = std::cpp20::ostream_iterator<char>;
using forward_iterator = std::forward_list<char>::iterator;
using bidirectional_iterator = std::list<char>::iterator;
using random_access_iterator = std::vector<char>::iterator;
using forward_iterator_const = std::forward_list<char>::const_iterator;
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
struct input_or_output_iter_value
{
    using type = std::iter_value_t<iterator_type>;
};

// ostream has std::iter_value_t void, but we want the output value type here.
template <typename value_t, typename ...ts>
struct input_or_output_iter_value<std::cpp20::ostream_iterator<value_t, ts...>>
{
    using type = value_t;
};

template <typename iterator_type>
using input_or_output_iter_value_t = typename input_or_output_iter_value<iterator_type>::type;

template <typename iterator_type>
struct test_sized_sentinel : public test_sentinel<input_or_output_iter_value_t<iterator_type>>
{
    using difference_type = typename iterator_type::difference_type;

    iterator_type pos;
};

template <typename iterator_t>
    requires std::random_access_iterator<iterator_t>
inline typename test_sized_sentinel<iterator_t>::difference_type
operator-(test_sized_sentinel<iterator_t> const & s,
          iterator_t const & i)
{
    return s.pos - i;
}

template <typename iterator_t>
    requires std::random_access_iterator<iterator_t>
inline typename test_sized_sentinel<iterator_t>::difference_type
operator-(iterator_t const & i,
          test_sized_sentinel<iterator_t> const & s)
{
    return i - s.pos;
}

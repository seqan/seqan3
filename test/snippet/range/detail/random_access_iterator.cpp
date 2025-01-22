// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/core/range/detail/random_access_iterator.hpp>

template <typename range_type>
class my_random_access_iterator :
    public seqan3::detail::random_access_iterator_base<range_type, seqan3::detail::random_access_iterator>
{
    //...
};

int main()
{}

// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/utility/char_operations/predicate.hpp>

int main()
{
    seqan3::is_in_interval<'A', 'G'>('C'); // returns true

    constexpr auto my_check = seqan3::is_in_interval<'A', 'G'>;
    my_check('H'); // returns false
}

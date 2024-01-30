// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/utility/char_operations/predicate.hpp>

int main()
{
    seqan3::is_char<'C'>('C'); // returns true

    constexpr auto my_check = seqan3::is_char<'C'>;
    my_check('c'); // returns false, because case is different
}

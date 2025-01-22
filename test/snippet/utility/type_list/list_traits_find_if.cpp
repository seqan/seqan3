// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/utility/type_list/traits.hpp>

int main()
{
    using list_t = seqan3::type_list<int, float, double>;

    // None of the types in list_t is a pointer so find_if returns -1. However, int and bool are both integral,
    // so find_if returns 0 for the first occurrence.
    static_assert(seqan3::list_traits::find_if<std::is_pointer, list_t> == -1);
    static_assert(seqan3::list_traits::find_if<std::is_integral, list_t> == 0);
}

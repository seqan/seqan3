// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/utility/type_pack/traits.hpp>

int main()
{
    // None of the types in t is a pointer so find_if returns -1. However, int and bool are both integral,
    // so find_if returns 0 for the first occurrence.
    static_assert(seqan3::pack_traits::find_if<std::is_pointer, int, float, double> == -1);
    static_assert(seqan3::pack_traits::find_if<std::is_integral, int, float, double> == 0);
}

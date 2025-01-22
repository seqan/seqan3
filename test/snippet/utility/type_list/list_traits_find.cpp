// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/utility/type_list/traits.hpp>

int main()
{
    using list_t = seqan3::type_list<int, float, bool>;

    // Double is not in list_t so find returns -1. However, bool is in the type list so find will return 2.
    static_assert(seqan3::list_traits::find<double, list_t> == -1);
    static_assert(seqan3::list_traits::find<bool, list_t> == 2);
}

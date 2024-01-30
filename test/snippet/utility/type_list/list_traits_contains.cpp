// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/utility/type_list/traits.hpp>

int main()
{
    using list_t = seqan3::type_list<int, float, bool, int>;

    static_assert(seqan3::list_traits::contains<double, list_t> == false); // Is type double in type_list list_t?

    static_assert(seqan3::list_traits::contains<float, list_t> == true); // Is type float in type_list list_t?
}

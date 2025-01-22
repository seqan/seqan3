// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/utility/type_pack/traits.hpp>

int main()
{
    using split_t = seqan3::pack_traits::split_after<3, int, float, bool, double, char, int>;

    // Use ::first_type and ::second_type to access the type lists after being split.
    static_assert(std::same_as<seqan3::type_list<int, float, bool>, split_t::first_type>);
    static_assert(std::same_as<seqan3::type_list<double, char, int>, split_t::second_type>);
}

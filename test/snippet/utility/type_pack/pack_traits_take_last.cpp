// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/utility/type_pack/traits.hpp>

int main()
{
    // Take the last two types in the pack.
    static_assert(std::same_as<seqan3::type_list<bool, int>, seqan3::pack_traits::take_last<2, int, float, bool, int>>);
}

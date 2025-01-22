// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/utility/type_pack/traits.hpp>

int main()
{
    // Return the a type list of the pack without the first type.
    static_assert(
        std::same_as<seqan3::type_list<float, bool, int>, seqan3::pack_traits::drop_front<int, float, bool, int>>);
}

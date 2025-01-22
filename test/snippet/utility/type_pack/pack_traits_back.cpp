// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/utility/type_pack/traits.hpp>

int main()
{
    // Check if the last value is float.
    static_assert(std::same_as<float, seqan3::pack_traits::back<int, float, bool, int, float>>);
}

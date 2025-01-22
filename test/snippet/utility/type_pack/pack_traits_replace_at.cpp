// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/utility/type_pack/traits.hpp>

int main()
{
    // Replace the second element of the type pack with int.
    static_assert(std::same_as<seqan3::type_list<int, int, bool, double>,
                               seqan3::pack_traits::replace_at<int, 1, int, float, bool, double>>);
}

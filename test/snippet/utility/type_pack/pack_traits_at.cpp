// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/utility/type_pack/traits.hpp>

int main()
{
    // Look at the 2nd element.
    static_assert(std::same_as<float, seqan3::pack_traits::at<1, int, float, bool, double>>);
    // Look at the last element.
    static_assert(std::same_as<double, seqan3::pack_traits::at<-1, int, float, bool, double>>);
}

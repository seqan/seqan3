// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/utility/type_list/traits.hpp>

int main()
{
    using list_t = seqan3::type_list<int, float, bool, double>;

    // Look at the 2nd element.
    static_assert(std::same_as<float, seqan3::list_traits::at<1, list_t>>);
    // Look at the last element.
    static_assert(std::same_as<double, seqan3::list_traits::at<-1, list_t>>);
}

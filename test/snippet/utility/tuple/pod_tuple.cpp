// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/tuple/pod_tuple.hpp>

int main()
{
    seqan3::pod_tuple<int, float> tuple1{3, 4.7};
    static_assert(std::is_standard_layout_v<seqan3::pod_tuple<int, float>>);
    static_assert(seqan3::trivial<seqan3::pod_tuple<int, float>>);
    seqan3::debug_stream << std::get<int>(tuple1) << '\n'; // 3

    // template parameters are automatically deduced:
    seqan3::pod_tuple tuple2{17, 3.7f, 19l};
    seqan3::debug_stream << std::get<0>(tuple2) << '\n'; // 17

    auto [i, f, l] = tuple2;                                   // creates an int i with value 17, float f...
    seqan3::debug_stream << i << ',' << f << ',' << l << '\n'; // 17,3.7,19
}

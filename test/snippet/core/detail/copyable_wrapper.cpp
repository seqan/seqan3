// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/core/detail/copyable_wrapper.hpp>

int main()
{
    int outer{};
    // Might be used for non-copyable lambdas. In this example, the lambda would be copyable even without the wrapper.
    seqan3::detail::copyable_wrapper wrapper{[&outer](int const x)
                                             {
                                                 outer += x;
                                                 return outer;
                                             }};
    auto wrapper_2 = wrapper;                     // Would not work with non-copyable lambda.
    seqan3::debug_stream << wrapper(2) << '\n';   // 2
    seqan3::debug_stream << wrapper_2(4) << '\n'; // 6
}

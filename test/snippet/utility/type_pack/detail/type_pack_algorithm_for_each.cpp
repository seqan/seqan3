// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <string>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/type_pack/detail/type_pack_algorithm.hpp>

int main()
{
    auto fn = [](auto && a)
    {
        seqan3::debug_stream << a;
    };

    // prints each argument, i.e. "0, 1, 2, 3\n"
    seqan3::detail::for_each(fn, 0, ", ", 1.0, ", ", std::string{"2, 3"}, '\n');

    // is the same as explicitly writing
    fn(0);
    fn(", ");
    fn(1.0);
    fn(", ");
    fn(std::string{"2, 3"});
    fn('\n');
    return 0;
}

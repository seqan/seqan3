// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/container/small_string.hpp>

int main()
{
    constexpr seqan3::small_string sm{"hello"}; // construct from string literal at compile time!

    static_assert(sm[0] == 'h'); // This way I can also test it at compile-time

    seqan3::debug_stream << sm.size() << '\n'; // prints 5! (the null character is only stored internally)

    // conversion to a normal string:
    std::string sm_string{sm.str()};
    // access data directly with a pointer to the underlying zero-terminated array:
    char const * sm_cstr{sm.c_str()};

    seqan3::debug_stream << sm << sm_string << sm_cstr << '\n'; // prints "hellohellohello"
}

// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alphabet/structure/dot_bracket3.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::dot_bracket3 letter{'.'_db3};

    letter.assign_char('(');
    seqan3::debug_stream << letter << '\n'; // prints "("

    letter.assign_char('F');                // Unknown characters are implicitly converted to '.'.
    seqan3::debug_stream << letter << '\n'; // prints "."
}

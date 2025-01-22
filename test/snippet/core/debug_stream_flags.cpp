// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/core/debug_stream.hpp>

int main()
{
    uint8_t i = 71;
    seqan3::debug_stream << '\'' << i << "'\n";                          // prints '71' (because flag is set by default)
    seqan3::debug_stream.unsetf(seqan3::fmtflags2::small_int_as_number); // unsets the flag
    seqan3::debug_stream << '\'' << i << "'\n";                          // prints 'G'
    seqan3::debug_stream << seqan3::fmtflags2::small_int_as_number << '\'' << i << "'\n"; // prints '71' again
    // instead of formatting the stream "inline", one can also call .setf()
}

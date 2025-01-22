// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <sstream>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using namespace seqan3::literals;

    std::ostringstream o;
    seqan3::debug_stream.set_underlying_stream(o);

    seqan3::debug_stream << "ACGT"_dna5;

    o.flush();
    seqan3::debug_stream << o.str(); // prints the string stream's buffer: "ACGT"
}

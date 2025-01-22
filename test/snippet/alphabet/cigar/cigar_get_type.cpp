// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alphabet/cigar/cigar.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using seqan3::get;
    using namespace seqan3::literals;

    seqan3::cigar letter{10, 'M'_cigar_operation};

    // Note that this is equivalent to get<0>(letter)
    uint32_t size{get<uint32_t>(letter)};

    // Note that this is equivalent to get<1>(letter)
    seqan3::cigar::operation cigar_char{get<seqan3::cigar::operation>(letter)};

    seqan3::debug_stream << "Size is " << size << '\n';
    seqan3::debug_stream << "Cigar char is " << cigar_char << '\n'; // seqan3::debug_stream converts to char on the fly.
}

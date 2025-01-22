// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <sstream>
#include <vector>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/input.hpp>

auto input = R"(>TEST1
ACGT
>Test2
AGGCTGA
>Test3
GGAGTATAATATATATATATATAT)";

int main()
{
    using seqan3::get;

    seqan3::sequence_file_input fin{std::istringstream{input},
                                    seqan3::format_fasta{},
                                    seqan3::fields<seqan3::field::id, seqan3::field::seq, seqan3::field::qual>{}};

    for (auto & [id, seq, qual] : fin) // the order is now different, "id" comes first, because it was specified first
    {
        seqan3::debug_stream << "ID:  " << id << '\n';
        seqan3::debug_stream << "SEQ: " << seq << '\n';
        seqan3::debug_stream << "QUAL: " << qual << '\n';
    }
}

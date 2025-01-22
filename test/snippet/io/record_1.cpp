// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <sstream>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/record.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/output.hpp>

auto input = R"(> TEST1
ACGT
> Test2
AGGCTGA
> Test3
GGAGTATAATATATATATATATAT)";

int main()
{
    // specify custom field combination/order to file:
    seqan3::sequence_file_input fin{std::istringstream{input},
                                    seqan3::format_fasta{},
                                    seqan3::fields<seqan3::field::id, seqan3::field::seq>{}};

    auto record = fin.front(); // get current record, in this case the first

    auto & id = record.id();
    seqan3::debug_stream << id << '\n'; // TEST1
    auto & seq = record.sequence();
    seqan3::debug_stream << seq << '\n'; // ACGT
}

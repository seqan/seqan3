// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <sstream>
#include <string>

#include <seqan3/alphabet/container/concatenated_sequences.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/io/sequence_file/format_fasta.hpp>
#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/utility/views/zip.hpp>

using namespace seqan3::literals;

struct data_storage_t
{
    seqan3::concatenated_sequences<seqan3::dna4_vector> sequences{"ACGT"_dna4, "AAA"_dna4};
    seqan3::concatenated_sequences<std::string> ids{std::string{"ID1"}, std::string{"ID2"}};
};

int main()
{
    data_storage_t data_storage{};

    // ... in your file writing function:

    seqan3::sequence_file_output fout{std::ostringstream{}, seqan3::format_fasta{}};

    fout = seqan3::views::zip(data_storage.sequences, data_storage.ids);
}

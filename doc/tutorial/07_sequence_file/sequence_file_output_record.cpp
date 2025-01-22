// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <string>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/io/sequence_file/record.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::sequence_file_output fout{std::cout, seqan3::format_fasta{}};

    using types = seqan3::type_list<std::vector<seqan3::dna5>, std::string>;
    using fields = seqan3::fields<seqan3::field::seq, seqan3::field::id>;
    using sequence_record_type = seqan3::sequence_record<types, fields>;

    for (int i = 0; i < 5; ++i) // ...
    {
        std::string id{"test_id"};
        seqan3::dna5_vector sequence{"ACGT"_dna5};

        sequence_record_type record{std::move(sequence), std::move(id)};

        fout.push_back(record);
    }
}

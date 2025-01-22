// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <sstream>
#include <string>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/io/sam_file/output.hpp>
#include <seqan3/utility/type_list/type_list.hpp>

int main()
{
    seqan3::sam_file_output fout{std::ostringstream{}, seqan3::format_sam{}};

    std::string ref_id;
    std::string read_id;

    std::vector<seqan3::dna5> read;

    // ... e.g. compute and alignment

    using alignment_type =
        std::pair<std::vector<seqan3::gapped<seqan3::dna5>>, std::vector<seqan3::gapped<seqan3::dna5>>>;

    alignment_type dummy_alignment{}; // an empty dummy alignment

    using types = seqan3::type_list<std::vector<seqan3::dna5>, std::string, alignment_type>;
    using types_as_ids = seqan3::fields<seqan3::field::seq, seqan3::field::id, seqan3::field::alignment>;
    // the record type specifies the fields we want to write
    using record_type = seqan3::record<types, types_as_ids>;

    // initialize record
    record_type rec{read, ref_id, dummy_alignment};

    // Write the record
    fout.push_back(rec);

    // same as
    fout.push_back(record_type{read, ref_id, dummy_alignment});

    // as all our fields are empty so this would print an
}

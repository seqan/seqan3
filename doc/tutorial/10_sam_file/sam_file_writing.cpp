// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/test/snippet/create_temporary_snippet_file.hpp>
// std::filesystem::current_path() / "out.sam" will be deleted after the execution
seqan3::test::create_temporary_snippet_file example_sam{"out.sam", ""};

//![main]
#include <seqan3/io/sam_file/all.hpp>

int main()
{
    using namespace seqan3::literals;

    auto filename = std::filesystem::current_path() / "out.sam";

    seqan3::sam_file_output fout{filename};

    using types = seqan3::type_list<std::vector<seqan3::dna5>, std::string, std::vector<seqan3::cigar>>;
    using fields = seqan3::fields<seqan3::field::seq, seqan3::field::id, seqan3::field::cigar>;
    using sam_record_type = seqan3::sam_record<types, fields>;

    // write the following to the file
    // r001	0	*	0	0	4M2I2M2D	*	0	0	ACGTACGT	*
    sam_record_type record{};
    record.id() = "r001";
    record.sequence() = "ACGTACGT"_dna5;
    record.cigar_sequence() = {{4, 'M'_cigar_operation},
                               {2, 'I'_cigar_operation},
                               {2, 'M'_cigar_operation},
                               {2, 'D'_cigar_operation}};

    fout.push_back(record);
}
//![main]

// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/test/snippet/create_temporary_snippet_file.hpp>

// std::filesystem::current_path() / "my.sam" will be deleted after the execution
seqan3::test::create_temporary_snippet_file mapping_sam{"my.sam", ""};

//![solution]
#include <filesystem>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/io/sam_file/all.hpp>

int main()
{
    using namespace seqan3::literals;

    auto filename = std::filesystem::current_path() / "my.sam";

    std::vector<std::string> ids{"read1", "read2"};
    std::vector<std::vector<seqan3::dna4>> seqs{"ACGATCGACTAGCTACGATCAGCTAGCAG"_dna4,
                                                "AGAAAGAGCGAGGCTATTTTAGCGAGTTA"_dna4};

    seqan3::sam_file_output fout{filename};

    using types = seqan3::type_list<std::string &, std::vector<seqan3::dna4> &, seqan3::sam_flag>;
    using fields = seqan3::fields<seqan3::field::id, seqan3::field::seq, seqan3::field::flag>;
    using sam_record_type = seqan3::sam_record<types, fields>;

    for (size_t i = 0; i < ids.size(); ++i)
    {
        fout.push_back(sam_record_type{ids[i], seqs[i], seqan3::sam_flag::unmapped});
    }
}
//![solution]

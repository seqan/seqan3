// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/test/snippet/create_temporary_snippet_file.hpp>
seqan3::test::create_temporary_snippet_file my_fastq{"my.fastq",
                                                     R"//![fastq_file](
@seq1
CGATCGATC
+
IIIIIIIII
@seq2
AGCG
+
IIII
@seq3
AGCTAGCAGCGATCG
+
IIIIIHIIJJIIIII
@seq4
AGC
+
III
@seq5
AGCTAGCAGCGATCG
+
IIIIIHIIJJIIIII
)//![fastq_file]"}; // std::filesystem::current_path() / "my.fastq" will be deleted after the execution

// std::filesystem::current_path() / "output.fastq" will be deleted after the execution
seqan3::test::create_temporary_snippet_file output{"output.fastq", ""};

//![solution]
#include <algorithm>
#include <filesystem>
#include <ranges>

#include <seqan3/io/sequence_file/all.hpp>

int main()
{
    std::filesystem::path current_path = std::filesystem::current_path();

    seqan3::sequence_file_input fin{current_path / "my.fastq"};
    seqan3::sequence_file_output fout{current_path / "output.fastq"};

    auto length_filter = std::views::filter(
        [](auto & rec)
        {
            return std::ranges::size(rec.sequence()) >= 5;
        });

    fout = fin | length_filter;

    // This would also work: copy an input range into an output range
    std::ranges::move(fin | length_filter, fout.begin());
}
//![solution]

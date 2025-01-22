// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/test/snippet/create_temporary_snippet_file.hpp>
seqan3::test::create_temporary_snippet_file my_fastq{"my.fastq",
                                                     R"//![fastq_file](
@seq1
AGCTAGCAGCGATCG
+
IIIIIHIIIIIIIII
@seq2
CGATCGATC
+
IIIIIIIII
@seq3
AGCGATCGAGGAATATAT
+
IIIIHHGIIIIHHGIIIH
)//![fastq_file]"}; // std::filesystem::current_path() / "my.fastq" will be deleted after the execution

//![main]
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/utility/views/zip.hpp>

int main()
{
    // for simplicity we take the same file
    seqan3::sequence_file_input fin1{std::filesystem::current_path() / "my.fastq"};
    seqan3::sequence_file_input fin2{std::filesystem::current_path() / "my.fastq"};

    for (auto && [rec1, rec2] : seqan3::views::zip(fin1, fin2)) // && is important!
    {                                                           // because seqan3::views::zip returns temporaries
        if (rec1.id() != rec2.id())
            throw std::runtime_error("Your pairs don't match.");
    }
}
//![main]

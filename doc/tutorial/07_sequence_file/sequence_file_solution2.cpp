// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/test/snippet/create_temporary_snippet_file.hpp>
seqan3::test::create_temporary_snippet_file my_fasta{"my.fasta",
                                                     R"//![fasta_file](
>seq1
AGCT
>seq2
CGATCGA
)//![fasta_file]"}; // std::filesystem::current_path() / "my.fasta" will be deleted after the execution

//![solution]
#include <filesystem>
#include <ranges> // std::ranges::copy

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>

int main()
{
    std::filesystem::path current_path = std::filesystem::current_path();

    seqan3::sequence_file_input fin{current_path / "my.fasta"};

    using record_type = decltype(fin)::record_type;
    std::vector<record_type> records{};

    // You can use a for loop:
    for (auto & record : fin)
    {
        records.push_back(std::move(record));
    }

    // But you can also do this:
    std::ranges::copy(fin, std::back_inserter(records));

    seqan3::debug_stream << records << '\n';
}
//![solution]

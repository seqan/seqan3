// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/test/snippet/create_temporary_snippet_file.hpp>
seqan3::test::create_temporary_snippet_file my_fastq{"my.fastq", "\n"};

//![main]
#include <seqan3/io/sequence_file/all.hpp>

int main()
{
    seqan3::sequence_file_input fin{std::filesystem::current_path() / "my.fastq"};
    using record_type = typename decltype(fin)::record_type;

    // Because `fin` is a range, we can access the first element by dereferencing fin.begin()
    record_type rec = *fin.begin();
}
//![main]

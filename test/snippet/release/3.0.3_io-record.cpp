// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

//![main]
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/input.hpp>

int main()
{
    seqan3::sequence_file_input fin{"my.fastq"};

    for (auto && record : fin)
    {
        seqan3::debug_stream << "id: " << record.id() << '\n';
        seqan3::debug_stream << "sequence: " << record.sequence() << '\n';
        seqan3::debug_stream << "base_qualities: " << record.base_qualities() << '\n';
    }
}
//![main]

#include <seqan3/test/snippet/create_temporary_snippet_file.hpp>
// std::filesystem::current_path() / "my.fastq" will be deleted after the execution
seqan3::test::create_temporary_snippet_file input_fastq{"my.fastq", "\n"};

// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

// this macro will rename the `main` function below as `sam_file_filename_construction`.
// we just want to show that the syntax does compile.
// We do this to show the usage of std::cin within a copyable snippet while being able to execute this snippet
// without std::cin.
#define main sam_file_filename_construction

//![main]
#include <filesystem>

#include <seqan3/io/sam_file/all.hpp>

int main()
{
    auto filename = std::filesystem::current_path() / "my.sam";

    seqan3::sam_file_input fin_from_filename{filename};

    seqan3::sam_file_input fin_from_stream{std::cin, seqan3::format_sam{}};

    return 0;
}
//![main]

// this makes the snippet an executable
#undef main
int main()
{} // do nothing

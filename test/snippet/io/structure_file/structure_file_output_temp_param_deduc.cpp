// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <filesystem>

#include <seqan3/io/structure_file/output.hpp>

int main()
{
    auto tmp_file = std::filesystem::temp_directory_path() / "my.dbn";

    seqan3::structure_file_output fout{tmp_file}; // Vienna format detected, std::ofstream opened for file

    std::filesystem::remove(tmp_file);
}

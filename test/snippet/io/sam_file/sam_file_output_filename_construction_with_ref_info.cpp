// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <filesystem>
#include <string>
#include <vector>

#include <seqan3/io/sam_file/output.hpp>

int main()
{
    auto tmp_file = std::filesystem::temp_directory_path() / "my.sam";

    std::vector<std::string> ref_ids{"ref1", "ref2"};
    std::vector<size_t> ref_lengths{1234, 5678};

    seqan3::sam_file_output fout{tmp_file, ref_ids, ref_lengths};
    std::filesystem::remove(tmp_file);
}

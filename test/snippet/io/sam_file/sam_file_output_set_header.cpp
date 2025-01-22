// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <sstream>
#include <string>
#include <vector>

#include <seqan3/io/sam_file/output.hpp>

int main()
{
    std::vector<std::string> ref_ids{"ref1", "ref2"};
    std::vector<size_t> ref_lengths{1234, 5678};

    // always give reference information if you want to have your header properly initialised
    seqan3::sam_file_output fout{std::ostringstream{}, ref_ids, ref_lengths, seqan3::format_sam{}};

    // add information to the header of the file.
    fout.header().comments.push_back("This is a comment");
}

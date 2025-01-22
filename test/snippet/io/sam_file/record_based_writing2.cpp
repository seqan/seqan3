// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <filesystem>
#include <sstream>
#include <tuple>

#include <seqan3/io/sam_file/output.hpp>

int main()
{
    // I only want to print the mapping position (field::ref_offset) and flag:
    seqan3::sam_file_output fout{std::ostringstream{},
                                 seqan3::format_sam{},
                                 seqan3::fields<seqan3::field::ref_offset, seqan3::field::flag>{}};

    unsigned mapping_pos{1300};
    seqan3::sam_flag flag{seqan3::sam_flag::none};

    // ...

    fout.emplace_back(mapping_pos, flag); // note that the order the arguments is now different, because
    // or:                                    you specified that REF_OFFSET should be first
    fout.push_back(std::tie(mapping_pos, flag));
}

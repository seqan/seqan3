// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <sstream>

#include <seqan3/io/sam_file/input.hpp>

auto input = R"(@HD	VN:1.6	SO:coordinate
r001	99	ref	7	30	8M2I4M1D3M	=	37	39	TTAGATAAAGGATACTG	*)";

int main()
{
    seqan3::sam_file_input fin{std::istringstream{input}, seqan3::format_sam{}};
    //                          ^ no need to specify the template arguments
}

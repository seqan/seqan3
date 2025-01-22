// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <sstream>

#include <seqan3/io/sam_file/all.hpp>

auto sam_file_raw = R"(First	0	*	0	0	*	*	0	0	ACGT	*
2nd	0	*	0	0	*	*	0	0	NATA	*
Third	0	*	0	0	*	*	0	0	GATA	*
)";

int main()
{
    // copying a file in one line:
    seqan3::sam_file_output{std::ostringstream{}, seqan3::format_sam{}} =
        seqan3::sam_file_input{std::istringstream{sam_file_raw}, seqan3::format_sam{}};

    // with seqan3::sam_file_output as a variable:
    seqan3::sam_file_output fout{std::ostringstream{}, seqan3::format_sam{}};
    seqan3::sam_file_input fin{std::istringstream{sam_file_raw}, seqan3::format_sam{}};
    fout = fin;

    // or in pipe notation:
    seqan3::sam_file_input{std::istringstream{sam_file_raw}, seqan3::format_sam{}}
        | seqan3::sam_file_output{std::ostringstream{}, seqan3::format_sam{}};
}

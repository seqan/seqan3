// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/core/platform.hpp>

//![snippet]
#include <ranges>
#include <sstream>

#include <seqan3/io/sam_file/all.hpp>

auto sam_file_raw = R"(@HD	VN:1.6	SO:coordinate	GO:none
@SQ	SN:ref	LN:45
r001	99	ref	7	30	*	=	37	39	TTAGATAAAGGATACTG	*
r003	0	ref	29	30	*	*	0	0	GCCTAAGCTAA	*	SA:Z:ref,29,-,6H5M,17,0;
r003	2064	ref	29	17	*	*	0	0	TAGGC	*	SA:Z:ref,9,+,5S6M,30,1;
r001	147	ref	237	30	*	=	7	-39	CAGCGGCAT	*	NM:i:1
)";

int main()
{
    auto input_file = seqan3::sam_file_input{std::istringstream{sam_file_raw}, seqan3::format_sam{}};
    input_file | std::views::take(3) // take only the first 3 records
        | seqan3::sam_file_output{std::ostringstream{}, seqan3::format_sam{}};
}
//![snippet]

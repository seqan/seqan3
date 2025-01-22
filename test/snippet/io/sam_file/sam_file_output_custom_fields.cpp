// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <filesystem>
#include <sstream>

#include <seqan3/io/sam_file/all.hpp>

auto sam_file_raw = R"(@HD	VN:1.6	SO:coordinate	GO:none
@SQ	SN:ref	LN:45
r001	99	ref	7	30	8M2I4M1D3M	=	37	39	TTAGATAAAGGATACTG	*
r003	0	ref	29	30	5S6M	*	0	0	GCCTAAGCTAA	*	SA:Z:ref,29,-,6H5M,17,0;
r003	2064	ref	29	17	6H5M	*	0	0	TAGGC	*	SA:Z:ref,9,+,5S6M,30,1;
r001	147	ref	237	30	9M	=	7	-39	CAGCGGCAT	*	NM:i:1
)";

int main()
{
    // fin uses custom fields, fout uses the default fields.
    seqan3::sam_file_input fin{std::istringstream{sam_file_raw},
                               seqan3::format_sam{},
                               seqan3::fields<seqan3::field::ref_offset, seqan3::field::flag>{}};
    // output doesn't have to match the configuration of the input
    seqan3::sam_file_output fout{std::ostringstream{}, seqan3::format_sam{}};

    for (auto & r : fin)
        fout.push_back(r); // copy all the records.
}

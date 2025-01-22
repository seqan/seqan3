// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <filesystem>
#include <sstream>

#include <seqan3/io/sam_file/input.hpp>

auto sam_file_raw = R"(@HD	VN:1.6	SO:coordinate	GO:none
@SQ	SN:ref	LN:45
r001	99	ref	7	30	8M2I4M1D3M	=	37	39	TTAGATAAAGGATACTG	*
r003	0	ref	29	30	5S6M	*	0	0	GCCTAAGCTAA	*	SA:Z:ref,29,-,6H5M,17,0;
r003	2064	ref	29	17	6H5M	*	0	0	TAGGC	*	SA:Z:ref,9,+,5S6M,30,1;
r001	147	ref	237	30	9M	=	7	-39	CAGCGGCAT	*	NM:i:1
)";

int main()
{
    // Create the temporary file.
    auto tmp_file = std::filesystem::temp_directory_path() / "my.sam";
    std::ofstream tmp_stream{tmp_file};
    tmp_stream << sam_file_raw;
    tmp_stream.close();

    seqan3::sam_file_input fin{tmp_file}; // SAM format assumed, regular std::ifstream taken as stream
    std::filesystem::remove(tmp_file);
}

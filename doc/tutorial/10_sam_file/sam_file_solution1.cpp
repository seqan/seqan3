// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/test/snippet/create_temporary_snippet_file.hpp>
seqan3::test::create_temporary_snippet_file my_sam{"my.sam",
                                                   R"//![sam_file](
@HD	VN:1.6	SO:coordinate
@SQ	SN:ref	LN:45
r001	99	ref	7	30	8M2I4M1D3M	=	37	39	TTAGATAAAGGATACTG	*
r003	0	ref	9	30	5S6M	*	0	0	GCCTAAGCTAA	*	SA:Z:ref,29,-,6H5M,17,0;
r004	0	ref	16	30	6M14N5M	*	0	0	ATAGCTTCAGC	*
r003	2064	ref	29	17	5M	*	0	0	TAGGC	*	SA:Z:ref,9,+,5S6M,30,1;
r001	147	ref	37	30	9M	=	7	-39	CAGCGGCAT	*	NM:i:1
)//![sam_file]"}; // std::filesystem::current_path() / "my.sam" will be deleted after the execution

//![solution]
#include <filesystem>
#include <ranges>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sam_file/all.hpp>

int main()
{
    auto filename = std::filesystem::current_path() / "my.sam";

    seqan3::sam_file_input fin{filename};

    double sum{};
    size_t count{};

    std::ranges::for_each(fin,
                          [&sum, &count](auto & record)
                          {
                              sum += record.mapping_quality();
                              ++count;
                          });

    seqan3::debug_stream << "Average: " << (sum / count) << '\n';
}
//![solution]

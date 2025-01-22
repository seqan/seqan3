// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/test/snippet/create_temporary_snippet_file.hpp>
seqan3::test::create_temporary_snippet_file example_sam{"example.sam",
                                                        R"//![sam_file](
@HD	VN:1.6	SO:coordinate
@SQ	SN:ref	LN:45
r001	99	ref	7	30	8M2I4M1D3M	=	37	39	TTAGATAAAGGATACTG	*
r003	0	ref	9	30	5S6M	*	0	0	GCCTAAGCTAA	*
r004	0	ref	16	30	6M14N5M	*	0	0	ATAGCTTCAGC	*
r003	2064	ref	29	17	5M	*	0	0	TAGGC	*
r001	147	ref	37	30	9M	=	7	-39	CAGCGGCAT	*	NM:i:1
)//![sam_file]"}; // std::filesystem::current_path() / "example.sam" will be deleted after the execution

//![main]
#include <seqan3/io/sam_file/all.hpp>

int main()
{
    auto filename = std::filesystem::current_path() / "example.sam";

    seqan3::sam_file_input fin{filename};

    // for (seqan3::sam_record record : fin) // this will copy the record
    for (auto && record : fin) // this will pass the record by reference (no copy)
    {
        seqan3::debug_stream << record.id() << '\n';
        seqan3::debug_stream << record.sequence() << '\n';
        seqan3::debug_stream << record.flag() << '\n' << '\n';
    }
}
//![main]

inline void sam_file_record_copy_for()
{
    seqan3::sam_file_input fin{""};

    // if this does not work, please update the comment above
    for (seqan3::sam_record record : fin)
        ;
}

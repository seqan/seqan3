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
    using namespace seqan3::literals;

    auto filename = std::filesystem::current_path() / "example.sam";

    std::vector<std::string> ref_ids{"ref"}; // list of one reference name
    std::vector<seqan3::dna5_vector> ref_sequences{"AGAGTTCGAGATCGAGGACTAGCGACGAGGCAGCGAGCGATCGAT"_dna5};

    seqan3::sam_file_input fin{filename, ref_ids, ref_sequences};

    for (auto & record : fin)
        seqan3::debug_stream << record.alignment() << '\n'; // Now you can print the whole alignment!
}
//![main]

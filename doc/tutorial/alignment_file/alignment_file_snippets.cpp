#include <fstream>
#include <seqan3/std/filesystem>

//![include_aligned_sequence]
#include <seqan3/alignment/aligned_sequence/aligned_sequence_concept.hpp>
//![include_aligned_sequence]

struct write_file_dummy_struct
{
    write_file_dummy_struct()
    {

auto file_raw = R"//![sam_file](
@HD	VN:1.6	SO:coordinate
@SQ	SN:ref	LN:45
r001	99	ref	7	30	8M2I4M1D3M	=	37	39	TTAGATAAAGGATACTG	*
r003	0	ref	9	30	5S6M	*	0	0	GCCTAAGCTAA	*
r004	0	ref	16	30	6M14N5M	*	0	0	ATAGCTTCAGC	*
r003	2064	ref	29	17	5M	*	0	0	TAGGC	*
r001	147	ref	37	30	9M	=	7	-39	CAGCGGCAT	*	NM:i:1
)//![sam_file]";

        std::ofstream file{std::filesystem::temp_directory_path()/"example.sam"};
        std::string str{file_raw};
        file << str.substr(1); // skip first newline
    }
};

write_file_dummy_struct go{};

//![main]
#include <seqan3/io/alignment_file/all.hpp>

using namespace seqan3;

int main()
{
//![main]

{
//![writing]
    auto filename = std::filesystem::temp_directory_path()/"out.sam";

    alignment_file_output fout{filename, fields<field::FLAG, field::MAPQ>{}};

    size_t mymapq{60};
    uint8_t flag{0};

    // ...

    fout.emplace_back(flag, mymapq);
    // or:
    fout.push_back(std::tie(flag, mymapq));
//![writing]
}

{
//![file_extensions]
    debug_stream << format_sam::file_extensions << std::endl; // prints [fastq,fq]
    format_sam::file_extensions.push_back("sm");
//![file_extensions]
}

{
/*
//![filename_construction]
    alignment_file_input fin_from_filename{"/tmp/my.sam"};

    alignment_file_input fin_from_stream{std::cin, format_sam{}};
//![filename_construction]
*/
}


{
//![read_custom_fields]
    auto filename = std::filesystem::temp_directory_path()/"example.sam";

    alignment_file_input fin{filename, fields<field::ID, field::SEQ, field::FLAG>{}};

    for (auto & [id, seq, flag /*order!*/] : fin)
    {
        debug_stream << id << std::endl;
        debug_stream << seq << std::endl;
        debug_stream << flag << std::endl;
    }
//![read_custom_fields]
}

{
//![alignments_without_ref]
    auto filename = std::filesystem::temp_directory_path()/"example.sam";

    alignment_file_input fin{filename, fields<field::ID, field::ALIGNMENT>{}};

    for (auto & [ id, alignment ] : fin)
    {
        debug_stream << id << ": " << get<1>(alignment) << std::endl;
    }
//![alignments_without_ref]
}

{
//![alignments_with_ref]
    auto filename = std::filesystem::temp_directory_path()/"example.sam";

    std::vector<std::string> ref_ids{"ref"}; // list of one reference name
    std::vector<dna5_vector> ref_sequences{"AGAGTTCGAGATCGAGGACTAGCGACGAGGCAGCGAGCGATCGAT"_dna5};

    alignment_file_input fin{filename, ref_ids, ref_sequences, fields<field::ALIGNMENT>{}};

    for (auto & [ alignment ] : fin)
    {
        debug_stream << alignment << std::endl; // Now you can print the whole alignment!
    }
//![alignments_with_ref]
}

//![main_end]
}
//![main_end]

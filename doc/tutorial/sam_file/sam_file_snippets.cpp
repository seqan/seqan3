#include <fstream>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/std/filesystem>

struct write_file_dummy_struct
{
    std::filesystem::path const file_path = std::filesystem::temp_directory_path()/"example.sam";

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

        std::ofstream file{file_path};
        std::string str{file_raw};
        file << str.substr(1); // skip first newline
    }

    ~write_file_dummy_struct()
    {
        std::error_code ec{};
        std::filesystem::remove(file_path, ec);

        if (ec)
            seqan3::debug_stream << "[WARNING] Could not delete " << file_path << ". " << ec.message() << '\n';
    }
};

write_file_dummy_struct go{};

//![main]
#include <seqan3/io/alignment_file/all.hpp>

int main()
{
//![main]

{
//![writing]
    auto filename = std::filesystem::temp_directory_path()/"out.sam";

    seqan3::alignment_file_output fout{filename, seqan3::fields<seqan3::field::flag, seqan3::field::mapq>{}};

    size_t mymapq{0};
    seqan3::sam_flag flag{seqan3::sam_flag::unmapped};

    // ...

    fout.emplace_back(flag, mymapq);
    // or:
    fout.push_back(std::tie(flag, mymapq));
//![writing]
}

{
//![file_extensions]
    seqan3::debug_stream << seqan3::format_sam::file_extensions << '\n'; // prints [fastq,fq]
    seqan3::format_sam::file_extensions.push_back("sm");
//![file_extensions]
}

{
/*
//![filename_construction]
    seqan3::alignment_file_input fin_from_filename{"/tmp/my.sam"};

    seqan3::alignment_file_input fin_from_stream{std::cin, seqan3::format_sam{}};
//![filename_construction]
*/
}


{
//![read_custom_fields]
    auto filename = std::filesystem::temp_directory_path()/"example.sam";

    seqan3::alignment_file_input fin{filename, seqan3::fields<seqan3::field::id,
                                                              seqan3::field::seq,
                                                              seqan3::field::flag>{}};

    for (auto & [id, seq, flag /*order!*/] : fin)
    {
        seqan3::debug_stream << id << '\n';
        seqan3::debug_stream << seq << '\n';
        seqan3::debug_stream << flag << '\n';
    }
//![read_custom_fields]
}

{
//![alignments_without_ref]
    auto filename = std::filesystem::temp_directory_path()/"example.sam";

    seqan3::alignment_file_input fin{filename, seqan3::fields<seqan3::field::id, seqan3::field::alignment>{}};

    for (auto & [ id, alignment ] : fin)
    {
        seqan3::debug_stream << id << ": " << std::get<1>(alignment) << '\n';
    }
//![alignments_without_ref]
}

{
//![alignments_with_ref]
    using seqan3::operator""_dna5;

    auto filename = std::filesystem::temp_directory_path()/"example.sam";

    std::vector<std::string> ref_ids{"ref"}; // list of one reference name
    std::vector<seqan3::dna5_vector> ref_sequences{"AGAGTTCGAGATCGAGGACTAGCGACGAGGCAGCGAGCGATCGAT"_dna5};

    seqan3::alignment_file_input fin{filename, ref_ids, ref_sequences, seqan3::fields<seqan3::field::alignment>{}};

    for (auto & [ alignment ] : fin)
    {
        seqan3::debug_stream << alignment << '\n'; // Now you can print the whole alignment!
    }
//![alignments_with_ref]
}
std::filesystem::remove(std::filesystem::temp_directory_path()/"out.sam");
//![main_end]
}
//![main_end]

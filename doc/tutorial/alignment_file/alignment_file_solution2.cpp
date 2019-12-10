#include <fstream>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/std/filesystem>

struct write_file_dummy_struct
{
    std::filesystem::path const tmp_path = std::filesystem::temp_directory_path();

    write_file_dummy_struct()
    {

auto fasta_file_raw = R"//![ref_file](
>chr1
ACAGCAGGCATCTATCGGCGGATCGATCAGGCAGGCAGCTACTGG
>chr2
ACAGCAGGCATCTATCGGCGGATCGATCAGGCAGGCAGCTACTGTAATGGCATCAAAATCGGCATG
)//![ref_file]";

auto file_raw = R"//![sam_file](
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:45
@SQ	SN:chr2	LN:66
r001	99	chr1	7	60	8M2I4M1D3M	=	37	39	TTAGATAAAGGATACTG	*
r003	0	chr1	9	60	5S6M	*	0	0	GCCTAAGCTAA	*
r004	0	chr2	16	60	6M14N5M	*	0	0	ATAGCTTCAGC	*
r003	2064	chr2	18	10	5M	*	0	0	TAGGC	*
)//![sam_file]";

        std::ofstream file{tmp_path/"mapping.sam"};
        std::string str{file_raw};
        file << str.substr(1); // skip first newline

        std::ofstream reffile{tmp_path/"reference.fasta"};
        std::string fasta_file_rawstr{fasta_file_raw};
        reffile << fasta_file_rawstr.substr(1); // skip first newline
    }

    ~write_file_dummy_struct()
    {
        std::error_code ec{};
        std::filesystem::path file_path{};

        file_path = tmp_path/"mapping.sam";
        std::filesystem::remove(file_path, ec);
        if (ec)
            seqan3::debug_stream << "[WARNING] Could not delete " << file_path << ". " << ec.message() << '\n';

        file_path = tmp_path/"reference.fasta";
        std::filesystem::remove(file_path, ec);
        if (ec)
            seqan3::debug_stream << "[WARNING] Could not delete " << file_path << ". " << ec.message() << '\n';

    }
};

write_file_dummy_struct go{};

//![solution]
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/alignment_file/all.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/alphabet/gap/gap.hpp>
#include <seqan3/std/filesystem>
#include <seqan3/std/ranges>

using seqan3::get;

struct my_traits : public seqan3::sequence_file_input_default_traits_dna
{
    //!\brief The container for sequences is now std::vector seqan3::concatenated_sequences.
    template <typename _sequence_container>
    using sequence_container_container = std::vector<_sequence_container>;
};

int main()
{
    std::filesystem::path tmp_dir = std::filesystem::temp_directory_path(); // get the temp directory

    // read in reference information
    seqan3::sequence_file_input<my_traits> reference_file{tmp_dir/"reference.fasta"};
    seqan3::concatenated_sequences<std::string> ref_ids = get<seqan3::field::id>(reference_file);
    std::vector<std::vector<seqan3::dna5>> ref_seqs = get<seqan3::field::seq>(reference_file);

    seqan3::alignment_file_input mapping_file{tmp_dir/"mapping.sam",
                                              ref_ids,
                                              ref_seqs,
                                              seqan3::fields<seqan3::field::id,
                                                             seqan3::field::ref_id,
                                                             seqan3::field::mapq,
                                                             seqan3::field::alignment>{}};

    auto mapq_filter = std::views::filter([] (auto & rec) { return seqan3::get<seqan3::field::mapq>(rec) >= 30; });

    for (auto & [id, ref_id, mapq, alignment] : mapping_file | mapq_filter)
    {
        auto & ref = std::get<0>(alignment);
        size_t sum_ref{};
        std::ranges::for_each(ref.begin(), ref.end(), [&sum_ref] (auto c) { if (c == seqan3::gap{}) ++sum_ref; });

        auto & read = std::get<0>(alignment);
        size_t sum_read{};
        std::ranges::for_each(read.begin(), read.end(), [&sum_read] (auto c) { if (c == seqan3::gap{}) ++sum_read; });

        seqan3::debug_stream << id << " mapped against " << ref_id << " with "
                             << sum_read << " gaps in the read sequence and "
                             << sum_ref  << " gaps in the reference sequence." << std::endl;
    }
}
//![solution]

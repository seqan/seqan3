#include <fstream>

#include <seqan3/std/filesystem>

struct write_file_dummy_struct
{
    write_file_dummy_struct()
    {

auto file_raw = R"//![sam_file](
@HD	VN:1.6	SO:coordinate
@SQ	SN:ref	LN:45
r001	99	ref	7	30	8M2I4M1D3M	=	37	39	TTAGATAAAGGATACTG	*
r003	0	ref	9	30	5S6M	*	0	0	GCCTAAGCTAA	*	SA:Z:ref,29,-,6H5M,17,0;
r004	0	ref	16	30	6M14N5M	*	0	0	ATAGCTTCAGC	*
r003	2064	ref	29	17	5M	*	0	0	TAGGC	*	SA:Z:ref,9,+,5S6M,30,1;
r001	147	ref	37	30	9M	=	7	-39	CAGCGGCAT	*	NM:i:1
)//![sam_file]";

        std::ofstream file{std::filesystem::temp_directory_path()/"my.sam"};
        std::string str{file_raw};
        file << str.substr(1); // skip first newline
    }
};

write_file_dummy_struct go{};

//![solution]
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/alignment_file/all.hpp>
#include <seqan3/range/views/get.hpp>
#include <seqan3/std/filesystem>
#include <seqan3/std/ranges>

int main()
{
    std::filesystem::path tmp_dir = std::filesystem::temp_directory_path(); // get the temp directory

    seqan3::alignment_file_input fin{tmp_dir/"my.sam", seqan3::fields<seqan3::field::mapq>{}};

    double sum{};
    size_t c{};

    std::ranges::for_each(fin.begin(), fin.end(), [&sum, &c] (auto & rec)
    {
        sum += seqan3::get<seqan3::field::mapq>(rec);
        ++c;
    });

    seqan3::debug_stream << "Average: " << (sum/c) << '\n';
}
//![solution]

#include <fstream>

#include <seqan3/std/filesystem>

struct write_file_dummy_struct
{
    write_file_dummy_struct()
    {

auto file_raw = R"//![fastq_file](
@seq1
AGCTAGCAGCGATCG
+
IIIIIHIIIIIIIII
@seq2
CGATCGATC
+
IIIIIIIII
@seq3
AGCGATCGAGGAATATAT
+
IIIIHHGIIIIHHGIIIH
)//![fastq_file]";

        std::ofstream file{std::filesystem::temp_directory_path()/"my.fastq"};
        std::string str{file_raw};
        file << str.substr(1); // skip first newline
    }
};

write_file_dummy_struct go{};

//![solution]
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/io/stream/debug_stream.hpp>
#include <seqan3/range/view/get.hpp>
#include <seqan3/std/filesystem>

using namespace seqan3;

int main()
{
    std::filesystem::path tmp_dir = std::filesystem::temp_directory_path(); // get the temp directory

    sequence_file_input fin{tmp_dir/"my.fastq", fields<field::ID, field::SEQ_QUAL>{}};

    for (auto & [id, seq_qual] : fin)
    {
        debug_stream << "ID:  "  << id << '\n';
        debug_stream << "SEQ: "  << (seq_qual | view::get<0>) << '\n';
        debug_stream << "QUAL: " << (seq_qual | view::get<1>) << '\n';
    }
}
//![solution]

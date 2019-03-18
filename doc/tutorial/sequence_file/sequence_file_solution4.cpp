#include <fstream>

#include <seqan3/std/filesystem>

struct write_file_dummy_struct
{
    write_file_dummy_struct()
    {

auto file_raw = R"//![fastq_file](
@seq1
CGATCGATC
+
IIIIII!!!
@seq2
AGCGATCGAGGAATATAT
+
IIIIHHGIII!!!!!!!!
@seq3
AGCTAGCAGCGATCG
+
IIIIIHII!!!!!!!
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
#include <seqan3/range/view/take_until.hpp>  // view::take_until
#include <seqan3/std/filesystem>

using namespace seqan3;

int main()
{
    std::filesystem::path tmp_dir = std::filesystem::temp_directory_path(); // get the temp directory

    sequence_file_input fin{tmp_dir/"my.fastq", fields<field::ID, field::SEQ_QUAL>{}};
    sequence_file_output fout{tmp_dir/"trimmed.fastq", fields<field::ID, field::SEQ_QUAL>{}};

    auto trimming = std::view::transform([] (auto & rec)
    {
        get<field::SEQ_QUAL>(rec) = get<field::SEQ_QUAL>(rec)
                                  | view::take_until([] (auto chr) { return chr.to_phred() <= 10; });
        return rec;
    });

    for (auto && rec : fin | trimming)
    {
        fout.push_back(rec);
    }
}
//![solution]

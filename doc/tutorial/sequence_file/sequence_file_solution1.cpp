#include <fstream>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/std/filesystem>

struct write_file_dummy_struct
{
    std::filesystem::path const file_path = std::filesystem::temp_directory_path()/"my.fastq";

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

//![solution]
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/std/filesystem>

int main()
{
    std::filesystem::path tmp_dir = std::filesystem::temp_directory_path(); // get the temp directory

    seqan3::sequence_file_input fin{tmp_dir/"my.fastq"};

    for (auto & rec : fin)
    {
        seqan3::debug_stream << "ID:  "  << seqan3::get<seqan3::field::id>(rec) << '\n';
        seqan3::debug_stream << "SEQ: "  << seqan3::get<seqan3::field::seq>(rec) << '\n';
        seqan3::debug_stream << "QUAL: " << seqan3::get<seqan3::field::qual>(rec) << '\n';
    }
}
//![solution]

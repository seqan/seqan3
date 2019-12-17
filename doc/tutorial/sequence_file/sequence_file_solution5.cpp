#include <fstream>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/std/filesystem>

struct write_file_dummy_struct
{
    std::filesystem::path const tmp_path = std::filesystem::temp_directory_path();

    write_file_dummy_struct()
    {

auto file_raw = R"//![fastq_file](
@seq1
CGATCGATC
+
IIIIIIIII
@seq2
AGCG
+
IIII
@seq3
AGCTAGCAGCGATCG
+
IIIIIHIIJJIIIII
@seq4
AGC
+
III
@seq5
AGCTAGCAGCGATCG
+
IIIIIHIIJJIIIII
)//![fastq_file]";

        std::ofstream file{tmp_path/"my.fastq"};
        std::string str{file_raw};
        file << str.substr(1); // skip first newline
    }

    ~write_file_dummy_struct()
    {
        std::error_code ec{};
        std::filesystem::path file_path{};

        file_path = tmp_path/"my.fastq";
        std::filesystem::remove(file_path, ec);
        if (ec)
            seqan3::debug_stream << "[WARNING] Could not delete " << file_path << ". " << ec.message() << '\n';

        file_path = tmp_path/"output.fastq";
        std::filesystem::remove(file_path, ec);
        if (ec)
            seqan3::debug_stream << "[WARNING] Could not delete " << file_path << ". " << ec.message() << '\n';

    }
};

write_file_dummy_struct go{};

//![solution]
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/range/views/persist.hpp>
#include <seqan3/std/filesystem>
#include <seqan3/std/ranges>

int main()
{
    std::filesystem::path tmp_dir = std::filesystem::temp_directory_path(); // get the temp directory

    seqan3::sequence_file_input fin{tmp_dir/"my.fastq"};
    seqan3::sequence_file_output fout{tmp_dir/"output.fastq"};

    auto length_filter = std::views::filter([] (auto & rec)
    {
        return std::ranges::size(seqan3::get<seqan3::field::seq>(rec)) >= 5;
    });

    fout = fin | length_filter;

    // This would also work:
    // fin | length_filter | fout;
}
//![solution]

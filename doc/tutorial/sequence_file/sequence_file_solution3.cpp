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

        std::ofstream file{std::filesystem::temp_directory_path()/"my.fastq"};
        std::string str{file_raw};
        file << str.substr(1); // skip first newline
    }
};

write_file_dummy_struct go{};

//![solution]
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/range/view/get.hpp>
#include <seqan3/std/filesystem>
#include <seqan3/std/ranges>

using namespace seqan3;

int main()
{
    std::filesystem::path tmp_dir = std::filesystem::temp_directory_path(); // get the temp directory

    sequence_file_input fin{tmp_dir/"my.fastq"};

    auto length_filter = std::view::filter([] (auto const & rec)
    {
        return std::ranges::size(get<field::SEQ>(rec)) >= 5;
    });

    // you can use a for loop

    // for (auto & rec : fin | length_filter | std::view::take(2))
    // {
    //     debug_stream << "ID: " << get<field::ID>(rec) << '\n';
    // }

    // But you can also do this :)
    std::vector<std::string> ids = std::move(fin | length_filter | std::view::take(2) | view::get<field::ID>);
    // Note that you need to know the type of id (std::string)
    // We use move to avoid copying

    debug_stream << ids << std::endl;
}
//![solution]

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
#include <seqan3/range/views/get.hpp>
#include <seqan3/range/views/move.hpp>
#include <seqan3/std/filesystem>
#include <seqan3/std/ranges>

int main()
{
    std::filesystem::path tmp_dir = std::filesystem::temp_directory_path(); // get the temp directory

    seqan3::sequence_file_input fin{tmp_dir/"my.fastq"};

    auto length_filter = std::views::filter([] (auto const & rec)
    {
        return std::ranges::size(seqan3::get<seqan3::field::seq>(rec)) >= 5;
    });

    // you can use a for loop

    // for (auto & rec : fin | length_filter | std::views::take(2))
    // {
    //     seqan3::debug_stream << "ID: " << seqan3::get<seqan3::field::id>(rec) << '\n';
    // }

    // But you can also do this to retrieve all IDs into a vector:
    std::vector<std::string> ids = fin
                                 | length_filter                                    // apply length filter
                                 | std::views::take(2)                              // take first two records
                                 | seqan3::views::get<seqan3::field::id>            // select only ID from record
                                 | seqan3::views::move                              // mark ID to be moved out of record
                                 | seqan3::views::to<std::vector<std::string>>;     // convert to container
    // Note that you need to know the type of id (std::string)

    seqan3::debug_stream << ids << '\n';
}
//![solution]

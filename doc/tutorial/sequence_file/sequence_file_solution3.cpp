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
#include <seqan3/std/filesystem>
#include <seqan3/std/ranges>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>

int main()
{
    std::filesystem::path tmp_dir = std::filesystem::temp_directory_path(); // get the temp directory

    seqan3::sequence_file_input fin{tmp_dir/"my.fastq"};

    auto length_filter = std::views::filter([] (auto const & rec)
    {
        return std::ranges::size(rec.sequence()) >= 5;
    });

    // Store all IDs into a vector:
    std::vector<std::string> ids{};
    for (auto & rec : fin | length_filter | std::views::take(2))
    {
        ids.push_back(std::move(rec.id()));
    }

    seqan3::debug_stream << ids << '\n';
}
//![solution]

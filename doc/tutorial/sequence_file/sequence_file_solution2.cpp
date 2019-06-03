#include <fstream>

#include <seqan3/std/filesystem>

struct write_file_dummy_struct
{
    write_file_dummy_struct()
    {

auto file_raw = R"//![fasta_file](
>seq1
AGCT
>seq2
CGATCGA
)//![fasta_file]";

        std::ofstream file{std::filesystem::temp_directory_path()/"my.fasta"};
        std::string str{file_raw};
        file << str.substr(1); // skip first newline
    }
};

write_file_dummy_struct go{};

//![solution]
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/std/filesystem>
#include <seqan3/std/ranges> // std::ranges::copy

using namespace seqan3;

int main()
{
    std::filesystem::path tmp_dir = std::filesystem::temp_directory_path(); // get the temp directory

    sequence_file_input fin{tmp_dir/"my.fasta"};

    using record_type = decltype(fin)::record_type;
    std::vector<record_type> records{};

    // You can use a for loop:

    // for (auto & rec : fin)
    // {
    //     records.push_back(std::move(rec));
    // }

    // But you can also do this:
    std::ranges::copy(fin, std::back_inserter(records));

    debug_stream << records << std::endl;
}
//![solution]

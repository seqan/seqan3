#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/std/filesystem>

int main()
{
    auto tmp_file = std::filesystem::temp_directory_path() / "my.fasta";

    seqan3::sequence_file_output fout{tmp_file}; // FastA format detected, std::ofstream opened for file

    std::filesystem::remove(tmp_file);
}

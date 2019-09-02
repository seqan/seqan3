#include <seqan3/io/alignment_file/output.hpp>
#include <seqan3/std/filesystem>

int main()
{
    auto tmp_file = std::filesystem::temp_directory_path() / "my.sam";

    seqan3::alignment_file_output fout{tmp_file}; // SAM format detected, std::ofstream opened for file

    std::filesystem::remove(tmp_file);
}

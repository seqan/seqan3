#include <seqan3/std/filesystem>

#include <seqan3/io/sam_file/output.hpp>

int main()
{
    auto tmp_file = std::filesystem::temp_directory_path() / "my.sam";

    seqan3::sam_file_output fout{tmp_file}; // SAM format detected, std::ofstream opened for file

    std::filesystem::remove(tmp_file);
}

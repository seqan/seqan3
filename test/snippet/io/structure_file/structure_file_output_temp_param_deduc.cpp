#include <seqan3/io/structure_file/output.hpp>
#include <seqan3/std/filesystem>

int main()
{
    auto tmp_file = std::filesystem::temp_directory_path() / "my.dbn";

    seqan3::structure_file_output fout{tmp_file}; // Vienna format detected, std::ofstream opened for file

    std::filesystem::remove(tmp_file);
}

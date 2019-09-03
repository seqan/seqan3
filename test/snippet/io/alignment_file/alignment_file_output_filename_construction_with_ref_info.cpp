#include <string>
#include <vector>

#include <seqan3/io/alignment_file/output.hpp>
#include <seqan3/std/filesystem>

int main()
{
    auto tmp_file = std::filesystem::temp_directory_path() / "my.sam";

    std::vector<std::string> ref_ids{"ref1", "ref2"};
    std::vector<size_t>      ref_lengths{1234, 5678};

    seqan3::alignment_file_output fout{tmp_file, ref_ids, ref_lengths};
    std::filesystem::remove(tmp_file);
}

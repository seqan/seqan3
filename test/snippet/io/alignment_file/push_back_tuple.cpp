//! [all]
#include <seqan3/io/alignment_file/all.hpp>

using namespace seqan3;

int main()
{
    alignment_file_output fout{std::filesystem::temp_directory_path()/"my.sam"};

    auto it = fout.begin();

    std::string id;
    dna5_vector seq;

    // ...

    fout.push_back(std::tie(seq, id));
}
//! [all]

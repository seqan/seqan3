//! [all]
#include <seqan3/io/alignment_file/all.hpp>

using namespace seqan3;

int main()
{
    alignment_file_output fout{filesystem::temp_directory_path()/"my.sam"};

    auto file_it = fout.begin();

    std::string id;
    dna5_vector seq;

    // ...

    // assign to file iterator
    *file_it = std::tie(seq, id);
    // is the same as:
    fout.push_back(std::tie(seq, id));
}
//! [all]

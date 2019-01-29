//! [all]
#include <seqan3/io/alignment_file/all.hpp>

using namespace seqan3;

int main()
{
    alignment_file_output fout{filesystem::temp_directory_path()/"my.sam"};

    auto it = fout.begin();

    record<type_list<uint32_t, std::string>, fields<field::MAPQ, field::ID>> r;

    // ...

    fout.push_back(r);
}
//! [all]

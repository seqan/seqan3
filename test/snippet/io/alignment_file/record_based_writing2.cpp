//! [all]
#include <seqan3/io/alignment_file/all.hpp>

using namespace seqan3;

int main()
{
    // I only want to print the mapping position (field::REF_OFFSET) and flag:
    alignment_file_output fout{std::filesystem::temp_directory_path()/"my.sam", fields<field::REF_OFFSET, field::FLAG>{}};

    unsigned mapping_pos{1300};
    unsigned flag{0};

    // ...

    fout.emplace_back(mapping_pos, flag);  // note that the order the arguments is now different, because
    // or:                                    you specified that REF_OFFSET should be first
    fout.push_back(std::tie(mapping_pos, flag));
}
//! [all]

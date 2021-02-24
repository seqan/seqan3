#include <sstream>

#include <seqan3/io/alignment_file/output.hpp>

int main()
{
    // I only want to print the mapping position (field::ref_offset) and flag:
    seqan3::alignment_file_output fout{std::ostringstream{},
                                       seqan3::format_sam{},
                                       seqan3::fields<seqan3::field::ref_offset, seqan3::field::flag>{}};

    unsigned mapping_pos{1300};
    seqan3::sam_flag flag{seqan3::sam_flag::none};

    // ...

    fout.emplace_back(mapping_pos, flag);  // note that the order the arguments is now different, because
    // or:                                    you specified that REF_OFFSET should be first
    fout.push_back(std::tie(mapping_pos, flag));
}

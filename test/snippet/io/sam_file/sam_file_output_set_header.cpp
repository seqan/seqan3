#include <sstream>
#include <string>
#include <vector>

#include <seqan3/io/sam_file/output.hpp>

int main()
{
    std::vector<std::string> ref_ids{"ref1", "ref2"};
    std::vector<size_t>      ref_lengths{1234, 5678};

    // always give reference information if you want to have your header properly initialised
    seqan3::sam_file_output fout{std::ostringstream{}, ref_ids, ref_lengths, seqan3::format_sam{}};

    // add information to the header of the file.
    fout.header().comments.push_back("This is a comment");
}

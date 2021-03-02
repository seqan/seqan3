#include <sstream>

#include <seqan3/io/sam_file/input.hpp>

auto input = R"(@HD	VN:1.6	SO:coordinate
r001	99	ref	7	30	8M2I4M1D3M	=	37	39	TTAGATAAAGGATACTG	*)";

int main()
{
    seqan3::sam_file_input fin{std::istringstream{input}, seqan3::format_sam{}};
    //                          ^ no need to specify the template arguments
}

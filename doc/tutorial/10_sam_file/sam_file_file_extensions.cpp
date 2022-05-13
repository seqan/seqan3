//![main]
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sam_file/all.hpp>

int main()
{
    seqan3::debug_stream << seqan3::format_sam::file_extensions << '\n'; // prints [sam]
    seqan3::format_sam::file_extensions.push_back("sm");
    seqan3::debug_stream << seqan3::format_sam::file_extensions << '\n'; // prints [sam,sm]
}
//![main]

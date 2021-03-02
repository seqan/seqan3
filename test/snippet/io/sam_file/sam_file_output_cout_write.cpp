#include <iostream>

#include <seqan3/io/sam_file/output.hpp>

int main()
{
    seqan3::sam_file_output fout{std::cout, seqan3::format_sam{}};
}

#include <iostream>

#include <seqan3/io/alignment_file/output.hpp>

int main()
{
    seqan3::alignment_file_output fout{std::cout, seqan3::format_sam{}};
}

#include <seqan3/io/sequence_file/output.hpp>

int main()
{
    seqan3::sequence_file_output fout{std::cout, seqan3::format_fasta{}};
    //                          ^ no need to specify the template arguments
}

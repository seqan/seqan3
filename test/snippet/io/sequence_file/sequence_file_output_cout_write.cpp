#include <iostream>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/io/sequence_file/output.hpp>

int main()
{
    using seqan3::operator""_dna5;

    seqan3::sequence_file_output fout{std::cout, seqan3::format_fasta{}};
    //                          ^ no need to specify the template arguments

    fout.emplace_back("ACGTN"_dna5, "example_id"); // default order for fasta: SEQ, ID
}

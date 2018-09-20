#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/io/alignment_file/all.hpp>
#include <iostream>

using namespace seqan3;
using namespace seqan3::literal;

int main()
{

{
//! [filename_construction]
alignment_file_output fout{"/tmp/my.sam"}; // SAM format detected, std::ofstream opened for file
//! [filename_construction]
}

{
//! [format_construction]
// no need to specify the template arguments <...> for format specialization:
alignment_file_output fout{"/tmp/my.sam", fields<field::MAPQ>{}};
//! [format_construction]
}

{
//! [write_range]
alignment_file_output fout{"/tmp/my.sam"};

std::vector<std::tuple<dna5_vector, std::string>> range
{
    { "ACGT"_dna5, "First" },
    { "NATA"_dna5, "2nd" },
    { "GATA"_dna5, "Third" }
}; // a range of "records"

fout = range; // will iterate over the records and write them
//! [write_range]
}

{
//! [set_header]
alignment_file_output fout{"/tmp/my.sam"};

// add information to the header of the file.
fout.header().comments.push_back("This is a comment");
//! [set_header]
}

}

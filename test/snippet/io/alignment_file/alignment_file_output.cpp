#include <iostream>

#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/io/alignment_file/all.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/range/view/persist.hpp>

using namespace seqan3;

int main()
{

// First create a /tmp/input.sam and /tmp/input.fasta
sequence_file_output fout{"/tmp/input.fasta"};
fout.emplace_back("ACGT"_dna4, "TEST1");
fout.emplace_back("AGGCTGA"_dna4, "Test2");
fout.emplace_back("ACTGA"_dna4, "Test2");

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

// \todo TODO: uncomment once seqan3::alignment_file_input is implemented.
#if 0
{
//! [custom_fields]
alignment_file_input  fin{"input.sam", fields<field::REF_OFFSET, field::FLAG>{}};
alignment_file_output fout{"output.sam"}; // doesn't have to match the configuration

for (auto & r : fin)
{
    fout.push_back(r); // copy all the records.
}
//! [custom_fields]
}
#endif

{
//! [input_range]
// file format conversion in one line:
alignment_file_output{"/tmp/output.sam"} = sequence_file_input{"/tmp/input.sam"};

// with alignment_file_output as a variable:
alignment_file_output fout{"/tmp/output.sam"};
sequence_file_input fin{"/tmp/input.fasta"};
fout = fin;

// or in pipe notation:
sequence_file_input{"/tmp/input.sam"} | alignment_file_output{"/tmp/output.sam"};
//! [input_range]
}

{
//! [io_pipeline]
sequence_file_input{"/tmp/input.sam"} | view::persist
                                      | ranges::view::take(5) // take only the first 5 records
                                      | alignment_file_output{"/tmp/output.sam"};
//! [io_pipeline]
}

{
//! [range]
alignment_file_output fout{"/tmp/my.sam"};

std::vector<std::tuple<dna5_vector, std::string>> range
{
    { "ACGT"_dna5, "First" },
    { "NATA"_dna5, "2nd" },
    { "GATA"_dna5, "Third" }
}; // a range of "records"

fout = range; // will iterate over the records and write them
//! [range]
}
}

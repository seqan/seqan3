//! [fileinput]
#include <string>                               // for std::string
#include <seqan3/io/sequence_file/input.hpp>    // for sequence_file_input
#include <seqan3/core/debug_stream.hpp>         // for debug_stream

int main ()
{   
    // Initialise a file input object with a FastA file.
    std::string filename = "seq.fasta";                 // filename: "seq.fasta"
    seqan3::sequence_file_input file_in{filename};   

    // Retrieve the sequences and ids.
    for (auto &[seq, id, qual] : file_in)
    {
        seqan3::debug_stream << "ID:     " << id << '\n';
        seqan3::debug_stream << "SEQ:    " << seq << '\n';
        seqan3::debug_stream << "Empty Qual." << qual << '\n';  // qual es empty for FastAfiles
    }

    return 0;
}
//! [fileinput]

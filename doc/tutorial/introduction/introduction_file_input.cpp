//! [fileinput_1]
#include <string>
#include <seqan3/core/debug_stream.hpp>         // for debug_stream
#include <seqan3/io/sequence_file/input.hpp>    // for sequence_file_input
//! [fileinput_1]

#include <filesystem>                           // to create temp directory 
#include <seqan3/alphabet/nucleotide/dna4.hpp>  // to create right datastructure in tmp file
#include <seqan3/io/sequence_file/output.hpp>   // to create tmp file

//! [fileinput_2]
int main ()
{
//! [fileinput_2]

    // Creating a tmp FastA file with sequences.
    // It's necessary, because some review processes in GitHub don't allow using
    // an alreadys existing file. So you can use the whole code from this file as
    // an example or the collection of snippets, shown in the introduction.

    auto tmp_dir = std::filesystem::temp_directory_path();
    std::string filename{tmp_dir/"seq.fasta"};
    {
        // Create a /tmp/my.fasta file.
        seqan3::sequence_file_output file_out{filename};

        using seqan3::operator""_dna4;  // for better readability of next code lines

        file_out.emplace_back("ACGTGATG"_dna4, std::string{"seq1"});
        file_out.emplace_back("AGTGATACT"_dna4, std::string{"seq2"});
    }

/*

// This is just to show the right source code in introduction (refering to collection of snippets).
//! [fileinput_3]
    // Initialise a file input object with a FastA file.
    std::string filename = "../input/seq.fasta";
//! [fileinput_3]
*/

//! [fileinput_4]
    seqan3::sequence_file_input file_in{filename};      // filename: "seq.fasta"

    // Retrieve the sequences and ids.
    for (auto &[seq, id, qual] : file_in)
    {
        seqan3::debug_stream << "ID:     " << id << '\n';
        seqan3::debug_stream << "SEQ:    " << seq << '\n';
        seqan3::debug_stream << "Empty Qual." << qual << '\n';  // qual es empty for FastAfiles
    }

    return 0;
}
//! [fileinput_4]
//! [read]
#include <string>                               // for std::string
#include <vector>                               // for std::vector

#include <seqan3/alphabet/nucleotide/dna5.hpp>  // for dna5 datastrucutre
#include <seqan3/argument_parser/all.hpp>       // for argument_parser
#include <seqan3/core/debug_stream.hpp>         // for debug_stream
#include <seqan3/io/sequence_file/input.hpp>    // for sequence_file_input

int main(int argc, char * argv[])
{
    // Receive the filename as program argument.
    std::string filename{};
    seqan3::argument_parser parser("My-first-program", argc, argv);
    parser.add_positional_option(filename, "The filename of the file to read.");

    try
    {
        parser.parse();
    }
    catch (seqan3::argument_parser_error const & ext)
    {
        seqan3::debug_stream << "[PARSER ERROR] " << ext.what() << '\n';
        return 0;
    }

    seqan3::debug_stream << "Reading file " << filename << "\n";

    // Create the vector to store sequences of type seqan3::dna5_vector.
    std::vector<seqan3::dna5_vector> sequences;

    // Iterate through the file and store the sequences.
    seqan3::sequence_file_input file_in{filename};
    for (auto & [ seq, id, qual ] : file_in)
    {
        sequences.push_back(seq);
    }

    seqan3::debug_stream << sequences << '\n';
    return 0;
}

//! [read]

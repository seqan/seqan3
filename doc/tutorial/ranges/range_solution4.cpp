#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/argument_parser/all.hpp>                       // include argument parser
#include <seqan3/range/container/bitcompressed_vector.hpp>      // include bitcompressed vector

using seqan3::operator""_dna4;

int main(int argc, char ** argv)
{
    seqan3::argument_parser myparser("Vector-implementations-comparison", argc, argv);
    size_t size{};
    bool use_bitvector{};
    myparser.add_positional_option(size, "Size of vector");
    myparser.add_flag(use_bitvector, 'b', "bitvector", "Use bitvector instead of vector");

    try
    {
         myparser.parse();
    }
    catch (seqan3::argument_parser_error const & ext)                     // catch user errors
    {
        seqan3::debug_stream << "[Error] " << ext.what() << "\n";
        return -1;
    }

    if (use_bitvector)
    {
        seqan3::bitcompressed_vector<seqan3::dna4> vector;
        vector.resize(size, 'A'_dna4);
        seqan3::debug_stream << "Allocated seqan3::bitcompressed_vector<seqan3::dna4> of size "
                             << vector.size() << '\n';
    }
    else
    {
        std::vector<seqan3::dna4> vector{size};
        seqan3::debug_stream << "Allocated std::vector<seqan3::dna4> of size " << vector.size() << '\n';
    }
}

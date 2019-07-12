#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/argument_parser/all.hpp>                       // include argument parser
#include <seqan3/range/container/bitcompressed_vector.hpp>      // include bitcompressed vector

using namespace seqan3;

int main(int argc, char ** argv)
{
    argument_parser myparser("Vector-implementations-comparison", argc, argv);
    size_t size{};
    bool use_bitvector{};
    myparser.add_positional_option(size, "Size of vector");
    myparser.add_flag(use_bitvector, 'b', "bitvector", "Use bitvector instead of vector");

    try
    {
         myparser.parse();
    }
    catch (parser_invalid_argument const & ext)                     // catch user errors
    {
        debug_stream << "[Error] " << ext.what() << "\n";
        return -1;
    }

    if (use_bitvector)
    {
        bitcompressed_vector<dna4> vector;
        vector.resize(size, 'A'_dna4);
        debug_stream << "Allocated bitcompressed_vector<dna4> of size " << vector.size() << std::endl;
    }
    else
    {
        std::vector<dna4> vector{size};
        debug_stream << "Allocated std::vector<dna4> of size " << vector.size() << std::endl;
    }
}

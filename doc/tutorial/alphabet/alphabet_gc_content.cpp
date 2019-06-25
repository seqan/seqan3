//! [exercise]
#include <array>     // std::array
#include <string>    // std::string
#include <vector>    // std::vector

#include <seqan3/alphabet/all.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/view/char_to.hpp>

using namespace seqan3;

int main (int argc, char * argv[])
{
    std::string input{};
    seqan3::argument_parser parser("GC-Content", argc, argv);
    parser.add_positional_option(input, "Specify an input sequence.");
    try
    {
        parser.parse();
    }
    catch (seqan3::parser_invalid_argument const & ext) // the input is invalid
    {
        debug_stream << "[PARSER ERROR] " << ext.what() << '\n';
        return 0;
    }

    // Convert the input to a dna5 sequence.
    std::vector<dna5> sequence{input | view::char_to<dna5>};

    // Initialise an array with count values for dna5 symbols.
    std::array<size_t, dna5::alphabet_size> count{}; // default initialised with zeroes

    // Increase the symbol count according to the sequence.
    for (dna5 symbol : sequence)
        ++count[symbol.to_rank()];

    // Calculate the GC content: (#G + #C) / (#A + #T + #G + #C).
    size_t gc = count['C'_dna5.to_rank()] + count['G'_dna5.to_rank()];
    size_t atgc = input.size() - count['N'_dna5.to_rank()];
    float gc_content = 1.0f * gc / atgc;

    debug_stream << "The GC content of " << sequence << " is " << 100 * gc_content << "%.\n";

    return 0;
}
//! [exercise]

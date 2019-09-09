#include <iostream>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/views/to_rank.hpp>

int main()
{
    using seqan3::operator""_dna5;

    // The alphabet normally needs to be converted to char explicitly:
    std::cout << seqan3::to_char('C'_dna5);                // prints 'C'

    // The debug_stream, on the other hand, does this automatically:
    seqan3::debug_stream << 'C'_dna5;                      // prints 'C'

    // The debug_stream can also print all types that model std::ranges::input_range:
    std::vector<seqan3::dna5> vec{"ACGT"_dna5};
    seqan3::debug_stream << vec;                           // prints "ACGT"

    // ranges of non-alphabets are printed comma-separated:
    seqan3::debug_stream << (vec | seqan3::views::to_rank); // prints "[0,1,2,3]"
}

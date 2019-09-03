#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/container/concatenated_sequences.hpp>

int main()
{
    using seqan3::operator""_dna4;
    
    seqan3::concatenated_sequences<seqan3::dna4_vector> foobar;
    foobar.insert(foobar.end(), 2, "ACGT"_dna4);
    seqan3::debug_stream << foobar[0] << '\n'; // "ACGT"
    seqan3::debug_stream << foobar[1] << '\n'; // "ACGT"
}

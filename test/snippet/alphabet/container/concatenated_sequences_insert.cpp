#include <seqan3/alphabet/container/concatenated_sequences.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::concatenated_sequences<seqan3::dna4_vector> foobar;
    foobar.insert(foobar.end(), "ACGT"_dna4);
    seqan3::debug_stream << foobar[0] << '\n'; // "ACGT"
}

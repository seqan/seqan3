#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/range/container/bitcompressed_vector.hpp>

int main()
{
    using namespace seqan3::literals;

    std::vector<seqan3::dna4>                  v0{"ACGT"_dna4}; // data occupies 4 bytes in memory
    seqan3::bitcompressed_vector<seqan3::dna4> v1{"ACGT"_dna4}; // data occupies 1 byte in memory
}

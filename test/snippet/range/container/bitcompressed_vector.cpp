#include <seqan3/range/container/bitcompressed_vector.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>

using namespace seqan3;
using namespace seqan3::literal;

int main()
{
//! [usage]
std::vector<dna4>          v0{"ACGT"_dna4}; // data occupies 4 bytes in memory
bitcompressed_vector<dna4> v1{"ACGT"_dna4}; // data occupies 1 byte in memory
//! [usage]
}

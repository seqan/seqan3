//! [usage]
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/view/kmer_hash.hpp>

using namespace seqan3;

int main()
{
    std::vector<dna4> text{"ACGTAGC"_dna4};
    std::vector<size_t> hashes = text | view::kmer_hash(3);
    debug_stream << hashes << '\n'; // [6,27,44,50,9]
}
//! [usage]

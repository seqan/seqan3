#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/view/kmer_hash.hpp>

int main()
{
    using seqan3::operator""_dna4;

    std::vector<seqan3::dna4> text{"ACGTAGC"_dna4};
    std::vector<size_t> hashes = text | seqan3::view::kmer_hash(3);
    seqan3::debug_stream << hashes << '\n'; // [6,27,44,50,9]
}

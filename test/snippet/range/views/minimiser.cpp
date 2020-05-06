#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/views/kmer_hash.hpp>
#include <seqan3/range/views/minimiser.hpp>

using seqan3::operator""_dna4;
using seqan3::operator""_shape;

int main()
{
    std::vector<seqan3::dna4> text{"ACGTAGC"_dna4};

    auto hashes = text | seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{3}});
    seqan3::debug_stream << hashes << '\n'; // [6,27,44,50,9]

    auto minimiser = hashes | seqan3::views::minimiser(4);
    seqan3::debug_stream << minimiser << '\n'; // [6,9]

    // kmer_hash with gaps, hashes: [2,7,8,14,1], minimiser: [2,1]
    seqan3::debug_stream << (text | seqan3::views::kmer_hash(0b101_shape) | seqan3::views::minimiser(4)) << '\n';
}

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/views/kmer_hash.hpp>

using namespace seqan3::literals;

int main()
{
    std::vector<seqan3::dna4> text{"ACGTAGC"_dna4};

    auto hashes = text | seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{3}});
    seqan3::debug_stream << hashes << '\n'; // [6,27,44,50,9]

    seqan3::debug_stream << (text | seqan3::views::kmer_hash(seqan3::ungapped{3})) << '\n'; // [6,27,44,50,9]

    seqan3::debug_stream << (text | seqan3::views::kmer_hash(0b101_shape)) << '\n'; // [2,7,8,14,1]
}

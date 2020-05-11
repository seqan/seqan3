#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/views/minimiser_hash.hpp>

using seqan3::operator""_dna4;
using seqan3::operator""_shape;

int main()
{
    std::vector<seqan3::dna4> text{"CCACGTCGACGGTT"_dna4};

    // Here a consecutive shape with size 4 (so the k-mer size is 4) and a window size of 8 is used. The seed is set
    // to 0, so lexicographical ordering is used for demonstration purposes.
    auto minimisers = text | seqan3::views::minimiser_hash(seqan3::shape{seqan3::ungapped{4}}, window_size{8}, seed{0});
    seqan3::debug_stream << minimisers << '\n'; // [27, 97, 26] representing the k-mers [ACGT, CGAC, ACGG]

    // Here a gapped shape with size 5 (and a k-mer size of 3) and a window size of 8 is used. The seed is set
    // to 0, so lexicographical ordering is used for demonstration purposes.
    auto minimisers2 = text | seqan3::views::minimiser_hash(0b10101_shape, window_size{8}, seed{0});
    seqan3::debug_stream << minimisers2 << '\n'; // [9, 18, 11] representing the k-mers [A.G.C, C.A.G, A.G.T]
}

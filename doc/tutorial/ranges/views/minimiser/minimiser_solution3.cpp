#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/views/kmer_hash.hpp>
#include <seqan3/range/views/minimiser_hash.hpp>

using seqan3::operator""_dna4;
using seqan3::operator""_shape;

int main()
{
    std::vector<seqan3::dna4> text{"CCACGTCGACGGTT"_dna4};

    // This would lead to an static assert error, because the shape.size() equals the window size. Therefore,
    // kmer_hash needs to be used.
    /*auto example_a = text | seqan3::views::minimiser_hash(seqan3::shape{seqan3::ungapped{4}},
                                                          seqan3::window_size{4}
                                                          seqan3::seed{0},
                                                          reverse = false);*/
    // results in: [81, 70, 27, 109, 182, 216, 97, 134, 26, 107, 175]
    // representing the k-mers [CCAC, CACG, ACGT, CGTC, GTCG, TCGA, CGAC, GACG, ACGG, CGGT, GGTT]
    auto example_a =  text | seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{4}});
    seqan3::debug_stream << example_a << '\n';

    auto example_b = text | seqan3::views::minimiser_hash(seqan3::shape{seqan3::ungapped{4}},
                                                          seqan3::window_size{8},
                                                          seqan3::seed{0},
                                                          reverse = false);
    // results in: [27, 97, 26] representing the k-mers [ACGT, CGAC, ACGG]
    seqan3::debug_stream << example_b << '\n';

    auto example_c = text | seqan3::views::minimiser_hash(0b10101_shape,
                                                          seqan3::window_size{8},
                                                          seqan3::seed{0},
                                                          reverse = false);
    // results in: [9, 18, 11] representing the k-mers [A.G.C, C.A.G, A.G.T]
    seqan3::debug_stream << example_c << '\n';
}

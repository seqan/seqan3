#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/core/simd/all.hpp>

using uint16x8_t = seqan3::simd_type_t<uint16_t, 8>;

int main()
{
    using seqan3::operator""_dna4;

    // Adds 7 sequences. The eighth will be defaulted.
    std::vector<seqan3::dna4_vector> batch;
    batch.push_back("ACGTACGTACGTACGATCG"_dna4);
    batch.push_back("AGTGAGCTACGGACTAGCTACGACT"_dna4);
    batch.push_back("GACTAGCACGAGCGAGATCG"_dna4);
    batch.push_back("GGATCGACGGACTAGC"_dna4);
    batch.push_back("ACGTACGACGGACGTACGAGCGAGCTACGAGC"_dna4);
    batch.push_back("ACGATCGACGACTAGCGAC"_dna4);
    batch.push_back("GTACGGATGGTAAACCGCACAT"_dna4);

    // Applies lazy transformation using `8` as a padding symbol if a sequence ends early.
    auto to_soa = batch | seqan3::views::to_simd<uint16x8_t>(8);

    size_t chunk_count = 0;
    for (auto && chunk : to_soa)
    {
        seqan3::debug_stream << "Chunk " << chunk_count++ << ":\n";
        for (auto & vec : chunk)
            seqan3::debug_stream << vec << '\n';

        seqan3::debug_stream << '\n';
    }
    return 0;
}

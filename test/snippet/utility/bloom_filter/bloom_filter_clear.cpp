#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/views/kmer_hash.hpp>
#include <seqan3/utility/bloom_filter/bloom_filter.hpp>

using seqan3::operator""_dna4;

int main()
{
    seqan3::bloom_filter bf{seqan3::bin_size{8192u},
                            seqan3::hash_function_count{2u}};

    auto const sequence1 = "ACTGACTGACTGATC"_dna4;
    auto const sequence2 = "GTGACTGACTGACTCG"_dna4;
    auto const sequence3 = "AAAAAAACGATCGACA"_dna4;
    auto hash_adaptor = seqan3::views::kmer_hash(seqan3::ungapped{5u});

    // Insert all 5-mers of sequence1
    for (auto && value : sequence1 | hash_adaptor)
        bf.emplace(value);

    // Insert all 5-mers of sequence3
    for (auto && value : sequence3 | hash_adaptor)
        bf.emplace(value);

    // Count all 5-mers of sequence2
    seqan3::debug_stream << bf.count(sequence2 | hash_adaptor) << '\n'; // 9

    // Clear the Bloom Filter
    bf.clear();

    // After clearing, no 5-mers are found
    seqan3::debug_stream << bf.count(sequence2 | hash_adaptor) << '\n'; // 0

}

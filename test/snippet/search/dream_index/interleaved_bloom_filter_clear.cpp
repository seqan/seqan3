#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <seqan3/search/views/kmer_hash.hpp>

using namespace seqan3::literals;

int main()
{
    seqan3::interleaved_bloom_filter ibf{seqan3::bin_count{8u},
                                         seqan3::bin_size{8192u},
                                         seqan3::hash_function_count{2u}};

    auto const sequence1 = "ACTGACTGACTGATC"_dna4;
    auto const sequence2 = "GTGACTGACTGACTCG"_dna4;
    auto const sequence3 = "AAAAAAACGATCGACA"_dna4;
    auto hash_adaptor = seqan3::views::kmer_hash(seqan3::ungapped{5u});

    // Insert all 5-mers of sequence1 into bin 0
    for (auto && value : sequence1 | hash_adaptor)
        ibf.emplace(value, seqan3::bin_index{0u});

    // Insert all 5-mers of sequence2 into bin 4
    for (auto && value : sequence2 | hash_adaptor)
        ibf.emplace(value, seqan3::bin_index{4u});

    // Insert all 5-mers of sequence3 into bin 7
    for (auto && value : sequence3 | hash_adaptor)
        ibf.emplace(value, seqan3::bin_index{7u});

    auto agent = ibf.counting_agent();

    // Count all 5-mers of sequence1 for all bins
    seqan3::debug_stream << agent.bulk_count(sequence1 | hash_adaptor) << '\n'; // [11,0,0,0,9,0,0,0]

    // Clear bin 0
    ibf.clear(seqan3::bin_index{0u});

    // After clearing, no 5-mers are found in bin 0
    seqan3::debug_stream << agent.bulk_count(sequence1 | hash_adaptor) << '\n'; // [0,0,0,0,9,0,0,0]

    // Search for specific values
    seqan3::debug_stream << agent.bulk_count(std::views::iota(0u, 1024u)) << '\n'; // [0,0,0,0,7,0,0,10]

    // Clear bin 4 and 7
    ibf.clear(std::vector{seqan3::bin_index{4u}, seqan3::bin_index{7u}});

    // After clearing, nothing is found
    seqan3::debug_stream << agent.bulk_count(std::views::iota(0u, 1024u)) << '\n'; // [0,0,0,0,0,0,0,0]
}

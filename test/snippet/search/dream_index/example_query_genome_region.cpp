#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/utility/views/chunk.hpp>
#include <seqan3/search/views/kmer_hash.hpp>

using namespace seqan3::literals;

int main()
{
    auto genome = "TTTTTTTTTTAAAAAAAAAATTTTTTTTTTGGGGGGGGGG"_dna4;

    seqan3::interleaved_bloom_filter ibf{seqan3::bin_count{4u}, seqan3::bin_size{8192u}}; // reserve 4 buckets

    // Fill buckets of the interleaved bloomm filter
    auto genome_buckets = genome | seqan3::views::chunk(10); // divide genome into buckets of size 10
    size_t bucket_idx{0};
    for (auto bucket : genome_buckets)
    {
        for (auto kmer : bucket | seqan3::views::kmer_hash(seqan3::ungapped{2})) // hash genome with k = 2
            ibf.emplace(kmer, seqan3::bin_index{bucket_idx});
        ++bucket_idx;
    }

    auto ibf_agent = ibf.counting_agent(); // the membership_agent enables efficient kemr queries

    auto query = "TTT"_dna4;
    auto query_kmers = query | seqan3::views::kmer_hash(seqan3::ungapped{2}); // hash query with k = 2

    seqan3::debug_stream << ibf_agent.bulk_count(query_kmers) << '\n'; // prints [2,0,2,0]
}

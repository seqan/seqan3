#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>

int main()
{
    // Construct an Interleaved Bloom Filter that contains 43 bins, each using 8192 bits, and 3 hash functions.
    seqan3::interleaved_bloom_filter ibf{seqan3::bin_count{43u},
                                         seqan3::bin_size{8192u},
                                         seqan3::hash_function_count{3}};

    // Construct an Interleaved Bloom Filter that contains 43 bins, each using 256 KiBits,
    // and the default of 2 hash functions.
    seqan3::interleaved_bloom_filter ibf2{seqan3::bin_count{43}, seqan3::bin_size{1ULL << 20}};
}

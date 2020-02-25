#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>

int main()
{
    // Construct an uncompressed Interleaved Bloom Filter.
    seqan3::interleaved_bloom_filter ibf{seqan3::bin_count{4u}, seqan3::bin_size{128u}, seqan3::hash_function_count{3}};

    // Fill `ibf` with data.

    // Construct an immutable, compressed Interleaved Bloom Filter.
    seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf2{ibf};
}

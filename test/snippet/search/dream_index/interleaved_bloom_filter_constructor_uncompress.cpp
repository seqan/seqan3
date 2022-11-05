#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>

int main()
{
    // Construct an uncompressed Interleaved Bloom Filter.
    seqan3::interleaved_bloom_filter ibf{seqan3::bin_count{9u}, seqan3::bin_size{128u}, seqan3::hash_function_count{3}};

    // Insert the values `126`, `712` and `237` into bins `0`, `3` and `8` of the Interleaved Bloom Filter.
    ibf.emplace(126, seqan3::bin_index{0u});
    ibf.emplace(712, seqan3::bin_index{3u});
    ibf.emplace(237, seqan3::bin_index{8u});

    // Construct an immutable, compressed Interleaved Bloom Filter.
    seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf_compressed{ibf};

    // Decompress the compressed Interleaved Bloom Filter.
    seqan3::interleaved_bloom_filter ibf_decompressed{ibf_compressed};
}

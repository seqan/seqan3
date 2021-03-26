#include <seqan3/utility/bloom_filter/bloom_filter.hpp>

int main()
{
    // Construct an uncompressed Bloom Filter.
    seqan3::bloom_filter bf{seqan3::bin_size{128u}, seqan3::hash_function_count{3}};

    // Fill `bf` with data.

    // Construct an immutable, compressed Bloom Filter.
    seqan3::bloom_filter<seqan3::data_layout::compressed> bf2{bf};
}

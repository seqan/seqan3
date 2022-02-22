#include <seqan3/utility/bloom_filter/bloom_filter.hpp>

int main()
{
    seqan3::bloom_filter bf{seqan3::bin_size{8192u}};

    // Insert the values `126`, `712` and `237` into the Bloom Filter.
    bf.emplace(126);
    bf.emplace(712);
    bf.emplace(237);
}

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>

int main()
{
    seqan3::interleaved_bloom_filter ibf{seqan3::bin_count{12u}, seqan3::bin_size{8192u}};
    ibf.emplace(126, seqan3::bin_index{0u});
    ibf.emplace(126, seqan3::bin_index{3u});
    ibf.emplace(126, seqan3::bin_index{9u});
    ibf.emplace(712, seqan3::bin_index{3u});
    ibf.emplace(237, seqan3::bin_index{9u});

    // The counting_vector must be at least as big as there are bins.
    seqan3::counting_vector<uint8_t> counts(12, 0);

    auto agent = ibf.membership_agent();

    counts += agent.bulk_contains(712);     // `counts` contains the number of occurrences of 712 in each bin.
    seqan3::debug_stream << counts << '\n'; // prints [0,0,0,1,0,0,0,0,0,0,0,0]

    counts += agent.bulk_contains(237);     // `counts` contains the number of occurrences of 712 and 237 in each bin.
    seqan3::debug_stream << counts << '\n'; // prints [0,0,0,1,0,0,0,0,0,1,0,0]

    counts += agent.bulk_contains(126); // `counts` contains the number of occurrences of 712, 237 and 126 in each bin.
    seqan3::debug_stream << counts << '\n'; // prints [1,0,0,2,0,0,0,0,0,2,0,0]

    counts += counts;                       // multiple counts can also be added together
    seqan3::debug_stream << counts << '\n'; // prints [2,0,0,4,0,0,0,0,0,4,0,0]
}

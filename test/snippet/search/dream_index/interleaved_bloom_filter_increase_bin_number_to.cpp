#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>

int main()
{
    seqan3::interleaved_bloom_filter ibf{seqan3::bin_count{12u}, seqan3::bin_size{8192u}};
    ibf.emplace(126, seqan3::bin_index{0u});
    ibf.emplace(712, seqan3::bin_index{3u});
    ibf.emplace(237, seqan3::bin_index{9u});

    ibf.increase_bin_number_to(seqan3::bin_count{18u});
    // Be sure to get the agent after `increase_bin_number_to` as it invalidates all agents!
    auto agent = ibf.membership_agent();

    // The content of the bins which were already present before the resize does not change
    seqan3::debug_stream << agent.bulk_contains(126) << '\n'; // prints [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    seqan3::debug_stream << agent.bulk_contains(712) << '\n'; // prints [0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    seqan3::debug_stream << agent.bulk_contains(237) << '\n'; // prints [0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0]
}

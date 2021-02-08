#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/views/kmer_hash.hpp>
#include <seqan3/search/dream_index/technical_binning_directory.hpp>

using seqan3::operator""_dna4;

int main()
{
    seqan3::ibf_config cfg{seqan3::bin_count{8u},
                           seqan3::bin_size{1ULL<<16},
                           seqan3::hash_function_count{2u}};

    std::vector<std::vector<seqan3::dna4>> technical_bins{"ACTGACTGACTGATC"_dna4,
                                                          "GTGACTGACTGACTCG"_dna4,
                                                          "AAAAAAACGATCGACA"_dna4};

    auto v = seqan3::views::kmer_hash(seqan3::ungapped{5u});

    seqan3::technical_binning_directory tbd2;

    seqan3::technical_binning_directory tbd{technical_bins, std::move(v), cfg};
    auto agent = tbd.counting_agent();
    seqan3::debug_stream << agent.count_query("ACTGACTGACTGATC"_dna4) << '\n'; // [11,9,0,0,0,0,0,0]

    seqan3::technical_binning_directory<seqan3::data_layout::compressed> ctbd{std::move(tbd)};
    auto agent2 = ctbd.counting_agent<uint16_t>(); // Also pass the type of the counter.
    seqan3::debug_stream << agent2.count_query("ACTGACTGACTGATC"_dna4) << '\n'; // [11,9,0,0,0,0,0,0]
}

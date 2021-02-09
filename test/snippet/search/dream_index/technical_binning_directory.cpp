#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/views/kmer_hash.hpp>
#include <seqan3/search/dream_index/technical_binning_directory.hpp>

#include <seqan3/test/tmp_filename.hpp>
#include <cereal/archives/binary.hpp>

using seqan3::operator""_dna4;

int main()
{
    seqan3::ibf_config cfg{seqan3::bin_count{8u},
                           seqan3::bin_size{1ULL<<16},
                           seqan3::hash_function_count{2u},
                        //    1u,
                           seqan3::ungapped{5u},
                           seqan3::window_size{19u},
                           seqan3::hash_variant::kmer};

    std::vector<std::vector<seqan3::dna4>> technical_bins{"ACTGACTGACTGATC"_dna4,
                                                          "GTGACTGACTGACTCG"_dna4,
                                                          "AAAAAAACGATCGACA"_dna4};

    seqan3::technical_binning_directory tbd{cfg};

    // auto v = seqan3::views::kmer_hash(seqan3::ungapped{5u});
    for (size_t i{0}; i < technical_bins.size(); ++i)
        tbd.emplace(technical_bins[i], seqan3::bin_index{i});

    auto agent = tbd.counting_agent();
    seqan3::debug_stream << agent.count_query("ACTGACTGACTGATC"_dna4) << '\n'; // [11,9,0,0,0,0,0,0]

    // seqan3::test::do_serialisation(tbd);

    using in_archive_t = cereal::BinaryInputArchive;
    using out_archive_t = cereal::BinaryOutputArchive;

    seqan3::test::tmp_filename filename{"cereal_test"};
    {
        std::ofstream os{filename.get_path(), std::ios::binary};
        out_archive_t oarchive{os};
        oarchive(tbd);
    }

    seqan3::technical_binning_directory tbd2;
    {
        std::ifstream is{filename.get_path(), std::ios::binary};
        in_archive_t iarchive{is};
        iarchive(tbd2);
    }

    auto agent2 = tbd2.counting_agent();
    seqan3::debug_stream << agent2.count_query("ACTGACTGACTGATC"_dna4) << '\n'; // [11,9,0,0,0,0,0,0]

    seqan3::technical_binning_directory<seqan3::data_layout::compressed> ctbd{std::move(tbd)};
    auto agent3 = ctbd.counting_agent<uint16_t>(); // Also pass the type of the counter.
    seqan3::debug_stream << agent3.count_query("ACTGACTGACTGATC"_dna4) << '\n'; // [11,9,0,0,0,0,0,0]
}

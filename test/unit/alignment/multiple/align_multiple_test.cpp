// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alignment/configuration/align_config_gap.hpp>
#include <seqan3/alignment/configuration/align_config_scoring_scheme.hpp>
#include <seqan3/alignment/scoring/aminoacid_scoring_scheme.hpp>
#include <seqan3/alignment/multiple/align_multiple.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/nucleotide/dna15.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/core/detail/debug_stream_alphabet.hpp>
#include <seqan3/range/views/to_char.hpp>
#include <seqan3/range/views/to.hpp>
#include <seqan3/range/views/deep.hpp>
#include <seqan3/test/expect_range_eq.hpp>

using seqan3::operator""_dna15;
using seqan3::operator""_dna4;
using seqan3::operator""_dna5;
using seqan3::operator""_rna5;
using seqan3::operator""_rna4;
using seqan3::operator""_aa27;

TEST(align_multiple_test, the_first_dna4_test)
{
    std::vector<seqan3::dna4_vector> const input{"AAAACCCGGG"_dna4, "AACCCGGG"_dna4, "AAAACGGG"_dna4};
    std::vector<std::string> const output{"AAAACCCGGG", "--AACCCGGG", "--AAAACGGG"};

    auto result = seqan3::align_multiple(input);
    auto result_strings = result | seqan3::views::deep{seqan3::views::to_char | seqan3::views::to<std::string>};

    EXPECT_RANGE_EQ(output, result_strings);
}

TEST(align_multiple_test, the_first_banded_test)
{
    std::vector<seqan3::dna4_vector> const input{"AAAACCCGGG"_dna4, "AACCCGGG"_dna4, "AAAACGGG"_dna4};
    std::vector<std::string> const output{"AAAACCCGGG", "--AACCCGGG", "--AAAACGGG"};

    auto cfg = seqan3::align_cfg::msa_default_configuration |
               seqan3::align_cfg::band_fixed_size{seqan3::align_cfg::lower_diagonal{-4},
                                                  seqan3::align_cfg::upper_diagonal{4}};

    auto result = seqan3::align_multiple(input, cfg);
    auto result_strings = result | seqan3::views::deep{seqan3::views::to_char | seqan3::views::to<std::string>};

    EXPECT_RANGE_EQ(output, result_strings);
}

TEST(align_multiple_test, the_first_aminoacid_test)
{
    // sequences taken from seqan/apps/seqan_tcoffee/tests/1aab.fa
    std::vector<seqan3::aa27_vector> const input
    {
        "KKDSNAPKRAMTSFMFFSSDFRSKHSDLSIVEMSKAAGAAWKELGPEERKVYEEMAEKDKERYKREM"_aa27,
        "KPKRPRSAYNIYVSESFQEAKDDSAQGKLKLVNEAWKNLSPEEKQAYIQLAKDDRIRYDNEMKSWEEQMAE"_aa27,
        "ADKPKRPLSAYMLWLNSARESIKRENPDFKVTEVAKKGGELWRGLKDKSEWEAKAATAKQNYIRALQEYERNGG"_aa27,
        "DPNKPKRAPSAFFVFMGEFREEFKQKNPKNKSVAAVGKAAGERWKSLSESEKAPYVAKANKLKGEYNKAIAAYNKGESA"_aa27
    };

    // alignment taken from seqan/apps/seqan_tcoffee/tests/1aab.fasta
    std::vector<std::string> const output
    {
        "KKDSNAPKRAMTSFMFFSSDFRSKHSDLS-----IVEMSKAAGAAWKELGPEERKVYEEMAEKDKERYKREM---------",
        "-----KPKRPRSAYNIYVSESFQEAKDDS-----AQGKLKLVNEAWKNLSPEEKQAYIQLAKDDRIRYDNEMKSWEEQMAE",
        "---ADKPKRPLSAYMLWLNSARESIKRENPDFK-VTEVAKKGGELWRGL--KDKSEWEAKAATAKQNYIRALQEYER-NGG",
        "--DPNKPKRAPSAFFVFMGEFREEFKQKNPKNKSVAAVGKAAGERWKSLSESEKAPYVAKANKLKGEYNKAIAAYNKGESA"
    };

    seqan3::configuration config = seqan3::align_cfg::scoring_scheme{seqan3::aminoacid_scoring_scheme{
                                                                        seqan3::aminoacid_similarity_matrix::BLOSUM62}};

    auto result = seqan3::align_multiple(input, config);
    auto result_strings = result | seqan3::views::deep{seqan3::views::to_char | seqan3::views::to<std::string>};

    EXPECT_RANGE_EQ(output, result_strings);
}

TEST(align_multiple_test, the_first_rna5_test)
{
    // alignment generated with seqan_tcoffee app: ./seqan_tcoffee -s data/rna5.fa -a rna -o data/out.fa
    std::vector<seqan3::rna5_vector> const input{"UUUNCCCGGG"_rna5, "UUCCCGGG"_rna5, "UUUNCGGG"_rna5};
    std::vector<std::string> const output{"UUUNCCCGGG","UU--CCCGGG","UU--UNCGGG"};

    auto result = seqan3::align_multiple(input);
    auto result_strings = result | seqan3::views::deep{seqan3::views::to_char | seqan3::views::to<std::string>};

    EXPECT_RANGE_EQ(output, result_strings);
}

TEST(align_multiple_test, the_first_gap_score_test)
{
    // alignment generated with seqan_tcoffee app: ./seqan_tcoffee -s data/dna1.fa -a dna -g -10 -e -2 -o data/out_g1.fa
    std::vector<seqan3::dna4_vector> const input{"ACGGTGG"_dna4, "ACCGTGCC"_dna4, "GCCGGTGCC"_dna4};
    std::vector<std::string> const output{"A-CGGTGG-", "A-CCGTGCC", "GCCGGTGCC"};

    constexpr seqan3::configuration cfg = seqan3::align_cfg::gap{seqan3::gap_scheme{seqan3::gap_score{-2},
                                                                                    seqan3::gap_open_score{-8}}};

    auto result = seqan3::align_multiple(input, cfg);
    auto result_strings = result | seqan3::views::deep{seqan3::views::to_char | seqan3::views::to<std::string>};

    EXPECT_RANGE_EQ(output, result_strings);
}

TEST(align_multiple_test, the_second_gap_score_test)
{
    // alignment generated with seqan_tcoffee app: ./seqan_tcoffee -s data/dna2.fa -a dna -g -2 -e -2 -o data/out_g2.fa
    std::vector<seqan3::dna5_vector> const input{"NNTGTNN"_dna5, "GGTNTNNGT"_dna5, "NGTNTGGG"_dna5};
    std::vector<std::string> const output{"NNTGTN-N-----", "G--GTNTNNGT--", "-N-GTNT--G-GG"};

    constexpr seqan3::configuration cfg = seqan3::align_cfg::gap{seqan3::gap_scheme{seqan3::gap_score{-2},
                                                                                    seqan3::gap_open_score{0}}};

    auto result = seqan3::align_multiple(input, cfg);
    auto result_strings = result | seqan3::views::deep{seqan3::views::to_char | seqan3::views::to<std::string>};

    EXPECT_RANGE_EQ(output, result_strings);
}

TEST(align_multiple_test, the_third_gap_score_test)
{
    // alignment generated with seqan_tcoffee app: ./seqan_tcoffee -s data/dna3.fa -a dna -g -16 -e -4 -o data/out_g3.fa
    std::vector<seqan3::dna15_vector> const input{"GGGTGGYTG"_dna15, "KTGTGGYTYTG"_dna15, "KTGTYYYTG"_dna15};
    std::vector<std::string> const output{"GGGTGG--YTG", "KTGTGGYTYTG", "KTGTYY--YTG"};

    constexpr seqan3::configuration cfg = seqan3::align_cfg::gap{seqan3::gap_scheme{seqan3::gap_score{-4},
                                                                                    seqan3::gap_open_score{-12}}};

    auto result = seqan3::align_multiple(input, cfg);
    auto result_strings = result | seqan3::views::deep{seqan3::views::to_char | seqan3::views::to<std::string>};

    EXPECT_RANGE_EQ(output, result_strings);
}

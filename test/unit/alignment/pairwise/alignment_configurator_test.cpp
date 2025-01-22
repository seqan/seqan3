// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <ranges>

#include <seqan3/alignment/configuration/align_config_score_type.hpp>
#include <seqan3/alignment/pairwise/alignment_configurator.hpp>
#include <seqan3/alignment/pairwise/detail/type_traits.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/test/expect_same_type.hpp>
#include <seqan3/utility/views/chunk.hpp>
#include <seqan3/utility/views/zip.hpp>

using seqan3::operator""_dna4;

auto setup()
{
    auto data = std::tuple{"ACGT"_dna4, "ACGT"_dna4};
    return std::views::single(std::move(data));
}

template <typename config_t>
auto run_test(config_t const & cfg)
{
    auto r = setup();
    auto [algorithm, complete_config] = seqan3::detail::alignment_configurator::configure<decltype(r)>(cfg);

    auto indexed_sequence_pairs = seqan3::views::zip(r, std::views::iota(0)) | seqan3::views::chunk(1);

    using complete_configuration_t = decltype(complete_config);
    using traits_t = seqan3::detail::alignment_configuration_traits<complete_configuration_t>;
    using alignment_result_t = typename traits_t::alignment_result_type;

    alignment_result_t align_result{};
    algorithm(*indexed_sequence_pairs.begin(),
              [&](auto && res) mutable
              {
                  align_result = std::forward<decltype(res)>(res);
              });

    return align_result;
}

TEST(alignment_configurator, configure_edit)
{
    EXPECT_EQ(run_test(seqan3::align_cfg::method_global{} | seqan3::align_cfg::edit_scheme).score(), 0);
}

TEST(alignment_configurator, configure_edit_end_position)
{
    EXPECT_EQ(run_test(seqan3::align_cfg::method_global{} | seqan3::align_cfg::edit_scheme
                       | seqan3::align_cfg::output_end_position{} | seqan3::align_cfg::output_score{})
                  .score(),
              0);
}

TEST(alignment_configurator, configure_edit_begin_position)
{
    EXPECT_EQ(run_test(seqan3::align_cfg::method_global{} | seqan3::align_cfg::edit_scheme
                       | seqan3::align_cfg::output_begin_position{} | seqan3::align_cfg::output_score{})
                  .score(),
              0);
}

TEST(alignment_configurator, configure_edit_trace)
{
    EXPECT_EQ(run_test(seqan3::align_cfg::method_global{} | seqan3::align_cfg::edit_scheme
                       | seqan3::align_cfg::output_alignment{} | seqan3::align_cfg::output_score{})
                  .score(),
              0);
}

TEST(alignment_configurator, configure_edit_semi)
{
    EXPECT_EQ(run_test(seqan3::align_cfg::method_global{seqan3::align_cfg::free_end_gaps_sequence1_leading{true},
                                                        seqan3::align_cfg::free_end_gaps_sequence2_leading{false},
                                                        seqan3::align_cfg::free_end_gaps_sequence1_trailing{true},
                                                        seqan3::align_cfg::free_end_gaps_sequence2_trailing{false}}
                       | seqan3::align_cfg::edit_scheme)
                  .score(),
              0);
}

TEST(alignment_configurator, configure_edit_banded)
{
    EXPECT_THROW((run_test(seqan3::align_cfg::method_global{} | seqan3::align_cfg::edit_scheme
                           | seqan3::align_cfg::band_fixed_size{seqan3::align_cfg::lower_diagonal{-1},
                                                                seqan3::align_cfg::upper_diagonal{1}})),
                 seqan3::invalid_alignment_configuration);
}

TEST(alignment_configurator, configure_edit_max_error)
{
    EXPECT_EQ(
        run_test(seqan3::align_cfg::method_global{} | seqan3::align_cfg::edit_scheme | seqan3::align_cfg::min_score{-3})
            .score(),
        0);
}

TEST(alignment_configurator, configure_affine_global)
{
    auto cfg =
        seqan3::align_cfg::method_global{}
        | seqan3::align_cfg::gap_cost_affine{seqan3::align_cfg::open_score{-10}, seqan3::align_cfg::extension_score{-1}}
        | seqan3::align_cfg::scoring_scheme{seqan3::nucleotide_scoring_scheme{}};

    EXPECT_EQ(run_test(cfg).score(), 0);
}

TEST(alignment_configurator, configure_affine_global_max_error)
{
    auto cfg =
        seqan3::align_cfg::method_global{}
        | seqan3::align_cfg::gap_cost_affine{seqan3::align_cfg::open_score{-10}, seqan3::align_cfg::extension_score{-1}}
        | seqan3::align_cfg::scoring_scheme{seqan3::nucleotide_scoring_scheme{}} | seqan3::align_cfg::min_score{-5};

    EXPECT_THROW(run_test(cfg), seqan3::invalid_alignment_configuration);
}

TEST(alignment_configurator, configure_affine_global_end_position)
{
    auto cfg =
        seqan3::align_cfg::method_global{}
        | seqan3::align_cfg::gap_cost_affine{seqan3::align_cfg::open_score{-10}, seqan3::align_cfg::extension_score{-1}}
        | seqan3::align_cfg::scoring_scheme{seqan3::nucleotide_scoring_scheme{}}
        | seqan3::align_cfg::output_end_position{} | seqan3::align_cfg::output_score{};

    EXPECT_EQ(run_test(cfg).score(), 0);
}

TEST(alignment_configurator, configure_affine_global_begin_position)
{
    auto cfg =
        seqan3::align_cfg::method_global{}
        | seqan3::align_cfg::gap_cost_affine{seqan3::align_cfg::open_score{-10}, seqan3::align_cfg::extension_score{-1}}
        | seqan3::align_cfg::scoring_scheme{seqan3::nucleotide_scoring_scheme{}}
        | seqan3::align_cfg::output_begin_position{} | seqan3::align_cfg::output_score{};

    EXPECT_EQ(run_test(cfg).score(), 0);
}

TEST(alignment_configurator, configure_affine_global_trace)
{
    auto cfg =
        seqan3::align_cfg::method_global{}
        | seqan3::align_cfg::gap_cost_affine{seqan3::align_cfg::open_score{-10}, seqan3::align_cfg::extension_score{-1}}
        | seqan3::align_cfg::scoring_scheme{seqan3::nucleotide_scoring_scheme{}} | seqan3::align_cfg::output_alignment{}
        | seqan3::align_cfg::output_score{};

    EXPECT_EQ(run_test(cfg).score(), 0);
}

TEST(alignment_configurator, configure_affine_global_banded)
{
    {
        auto cfg = seqan3::align_cfg::method_global{}
                 | seqan3::align_cfg::scoring_scheme{seqan3::nucleotide_scoring_scheme{}}
                 | seqan3::align_cfg::gap_cost_affine{seqan3::align_cfg::open_score{-10},
                                                      seqan3::align_cfg::extension_score{-1}}
                 | seqan3::align_cfg::band_fixed_size{seqan3::align_cfg::lower_diagonal{-1},
                                                      seqan3::align_cfg::upper_diagonal{1}};

        EXPECT_EQ(run_test(cfg).score(), 0);
    }

    { // invalid band
        auto cfg_base = seqan3::align_cfg::method_global{}
                      | seqan3::align_cfg::scoring_scheme{seqan3::nucleotide_scoring_scheme{}}
                      | seqan3::align_cfg::gap_cost_affine{seqan3::align_cfg::open_score{-10},
                                                           seqan3::align_cfg::extension_score{-1}};
        auto cfg_lower = cfg_base
                       | seqan3::align_cfg::band_fixed_size{seqan3::align_cfg::lower_diagonal{-10},
                                                            seqan3::align_cfg::upper_diagonal{-5}};
        auto cfg_upper = cfg_base
                       | seqan3::align_cfg::band_fixed_size{seqan3::align_cfg::lower_diagonal{5},
                                                            seqan3::align_cfg::upper_diagonal{6}};

        EXPECT_THROW(run_test(cfg_lower), seqan3::invalid_alignment_configuration);
        EXPECT_THROW(run_test(cfg_upper), seqan3::invalid_alignment_configuration);
    }
}

TEST(alignment_configurator, configure_affine_global_banded_with_alignment)
{
    auto cfg =
        seqan3::align_cfg::method_global{}
        | seqan3::align_cfg::gap_cost_affine{seqan3::align_cfg::open_score{-10}, seqan3::align_cfg::extension_score{-1}}
        | seqan3::align_cfg::scoring_scheme{seqan3::nucleotide_scoring_scheme{}}
        | seqan3::align_cfg::band_fixed_size{seqan3::align_cfg::lower_diagonal{-1},
                                             seqan3::align_cfg::upper_diagonal{1}};

    auto cfg_trace = cfg | seqan3::align_cfg::output_alignment{} | seqan3::align_cfg::output_score{};
    auto cfg_begin = cfg | seqan3::align_cfg::output_begin_position{} | seqan3::align_cfg::output_score{};
    auto cfg_end = cfg | seqan3::align_cfg::output_end_position{} | seqan3::align_cfg::output_score{};

    EXPECT_EQ(run_test(cfg_end).score(), 0);
    EXPECT_EQ(run_test(cfg_trace).score(), 0);
    EXPECT_EQ(run_test(cfg_begin).score(), 0);
}

TEST(alignment_configurator, configure_affine_global_semi)
{
    auto cfg = seqan3::align_cfg::method_global{seqan3::align_cfg::free_end_gaps_sequence1_leading{true},
                                                seqan3::align_cfg::free_end_gaps_sequence2_leading{true},
                                                seqan3::align_cfg::free_end_gaps_sequence1_trailing{true},
                                                seqan3::align_cfg::free_end_gaps_sequence2_trailing{true}}
             | seqan3::align_cfg::scoring_scheme{seqan3::nucleotide_scoring_scheme{}}
             | seqan3::align_cfg::gap_cost_affine{seqan3::align_cfg::open_score{-10},
                                                  seqan3::align_cfg::extension_score{-1}};

    EXPECT_EQ(run_test(cfg).score(), 0);
}

TEST(alignment_configurator, configure_affine_local)
{
    auto cfg =
        seqan3::align_cfg::method_local{}
        | seqan3::align_cfg::gap_cost_affine{seqan3::align_cfg::open_score{-10}, seqan3::align_cfg::extension_score{-1}}
        | seqan3::align_cfg::scoring_scheme{seqan3::nucleotide_scoring_scheme{}} | seqan3::align_cfg::output_score{};

    EXPECT_EQ(run_test(cfg).score(), 0);
}

TEST(alignment_configurator, configure_affine_local_end_positions)
{
    auto cfg =
        seqan3::align_cfg::method_local{}
        | seqan3::align_cfg::gap_cost_affine{seqan3::align_cfg::open_score{-10}, seqan3::align_cfg::extension_score{-1}}
        | seqan3::align_cfg::scoring_scheme{seqan3::nucleotide_scoring_scheme{}}
        | seqan3::align_cfg::output_end_position{} | seqan3::align_cfg::output_score{};

    EXPECT_EQ(run_test(cfg).score(), 0);
}

TEST(alignment_configurator, configure_affine_local_begin_positions)
{
    auto cfg =
        seqan3::align_cfg::method_local{}
        | seqan3::align_cfg::gap_cost_affine{seqan3::align_cfg::open_score{-10}, seqan3::align_cfg::extension_score{-1}}
        | seqan3::align_cfg::scoring_scheme{seqan3::nucleotide_scoring_scheme{}}
        | seqan3::align_cfg::output_begin_position{} | seqan3::align_cfg::output_score{};

    EXPECT_EQ(run_test(cfg).score(), 0);
}

TEST(alignment_configurator, configure_affine_local_alignment)
{
    auto cfg =
        seqan3::align_cfg::method_local{}
        | seqan3::align_cfg::gap_cost_affine{seqan3::align_cfg::open_score{-10}, seqan3::align_cfg::extension_score{-1}}
        | seqan3::align_cfg::scoring_scheme{seqan3::nucleotide_scoring_scheme{}} | seqan3::align_cfg::output_alignment{}
        | seqan3::align_cfg::output_score{};

    EXPECT_EQ(run_test(cfg).score(), 0);
}

TEST(alignment_configurator, configure_result_score_type)
{
    auto cfg = seqan3::align_cfg::method_global{} | seqan3::align_cfg::edit_scheme
             | seqan3::align_cfg::output_end_position{} | seqan3::align_cfg::output_score{}
             | seqan3::align_cfg::score_type<double>{};
    auto result = run_test(cfg);

    EXPECT_DOUBLE_EQ(result.score(), 0.0);
    EXPECT_EQ(result.sequence1_end_position(), 4u);
    EXPECT_EQ(result.sequence2_end_position(), 4u);
    EXPECT_SAME_TYPE(decltype(result.score()), double);
}

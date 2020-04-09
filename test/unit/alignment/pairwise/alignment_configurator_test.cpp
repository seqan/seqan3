// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <range/v3/view/single.hpp>

#include <seqan3/alignment/pairwise/alignment_configurator.hpp>
#include <seqan3/alignment/pairwise/detail/type_traits.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/range/views/chunk.hpp>
#include <seqan3/range/views/zip.hpp>
#include <seqan3/std/ranges>

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
    algorithm(*indexed_sequence_pairs.begin(), [&] (auto && res) mutable
    {
        align_result = std::forward<decltype(res)>(res);
    });

    return align_result;
}

TEST(alignment_configurator, configure_edit)
{
    EXPECT_EQ(run_test(seqan3::align_cfg::edit).score(), 0);
}

TEST(alignment_configurator, configure_edit_end_position)
{
    EXPECT_EQ(run_test(seqan3::align_cfg::edit | seqan3::align_cfg::result{seqan3::with_back_coordinate}).score(), 0);
}

TEST(alignment_configurator, configure_edit_begin_position)
{
    EXPECT_EQ(run_test(seqan3::align_cfg::edit | seqan3::align_cfg::result{seqan3::with_front_coordinate}).score(), 0);
}

TEST(alignment_configurator, configure_edit_trace)
{
    EXPECT_EQ(run_test(seqan3::align_cfg::edit | seqan3::align_cfg::result{seqan3::with_alignment}).score(), 0);
}

TEST(alignment_configurator, configure_edit_semi)
{
    EXPECT_EQ(run_test(seqan3::align_cfg::edit | seqan3::align_cfg::aligned_ends{seqan3::free_ends_first}).score(), 0);
}

TEST(alignment_configurator, configure_edit_banded)
{
    EXPECT_THROW((run_test(seqan3::align_cfg::edit |
                           seqan3::align_cfg::band{seqan3::static_band{seqan3::lower_bound{-1},
                                                                       seqan3::upper_bound{1}}})),
                 seqan3::invalid_alignment_configuration);
}

TEST(alignment_configurator, configure_edit_max_error)
{
    EXPECT_EQ(run_test(seqan3::align_cfg::edit | seqan3::align_cfg::max_error{3u}).score(), 0);
}

TEST(alignment_configurator, configure_affine_global)
{
    auto cfg = seqan3::align_cfg::mode{seqan3::global_alignment} |
               seqan3::align_cfg::gap{seqan3::gap_scheme{seqan3::gap_score{-1}, seqan3::gap_open_score{-10}}} |
               seqan3::align_cfg::scoring{seqan3::nucleotide_scoring_scheme{}};

    EXPECT_EQ(run_test(cfg).score(), 0);
}

TEST(alignment_configurator, configure_affine_global_max_error)
{
    auto cfg = seqan3::align_cfg::mode{seqan3::global_alignment} |
               seqan3::align_cfg::gap{seqan3::gap_scheme{seqan3::gap_score{-1}, seqan3::gap_open_score{-10}}} |
               seqan3::align_cfg::scoring{seqan3::nucleotide_scoring_scheme{}} |
               seqan3::align_cfg::max_error{5u};

    EXPECT_THROW(run_test(cfg), seqan3::invalid_alignment_configuration);
}

TEST(alignment_configurator, configure_affine_global_end_position)
{
    auto cfg = seqan3::align_cfg::mode{seqan3::global_alignment} |
               seqan3::align_cfg::gap{seqan3::gap_scheme{seqan3::gap_score{-1}, seqan3::gap_open_score{-10}}} |
               seqan3::align_cfg::scoring{seqan3::nucleotide_scoring_scheme{}} |
               seqan3::align_cfg::result{seqan3::with_back_coordinate};

    EXPECT_EQ(run_test(cfg).score(), 0);
}

TEST(alignment_configurator, configure_affine_global_begin_position)
{
    auto cfg = seqan3::align_cfg::mode{seqan3::global_alignment} |
               seqan3::align_cfg::gap{seqan3::gap_scheme{seqan3::gap_score{-1}, seqan3::gap_open_score{-10}}} |
               seqan3::align_cfg::scoring{seqan3::nucleotide_scoring_scheme{}} |
               seqan3::align_cfg::result{seqan3::with_front_coordinate};

    EXPECT_EQ(run_test(cfg).score(), 0);
}

TEST(alignment_configurator, configure_affine_global_trace)
{
    auto cfg = seqan3::align_cfg::mode{seqan3::global_alignment} |
               seqan3::align_cfg::gap{seqan3::gap_scheme{seqan3::gap_score{-1}, seqan3::gap_open_score{-10}}} |
               seqan3::align_cfg::scoring{seqan3::nucleotide_scoring_scheme{}} |
               seqan3::align_cfg::result{seqan3::with_alignment};

    EXPECT_EQ(run_test(cfg).score(), 0);
}

TEST(alignment_configurator, configure_affine_global_banded)
{
    {
        auto cfg = seqan3::align_cfg::mode{seqan3::global_alignment} |
                   seqan3::align_cfg::scoring{seqan3::nucleotide_scoring_scheme{}} |
                   seqan3::align_cfg::gap{seqan3::gap_scheme{seqan3::gap_score{-1}, seqan3::gap_open_score{-10}}} |
                   seqan3::align_cfg::band{seqan3::static_band{seqan3::lower_bound{-1}, seqan3::upper_bound{1}}};

        EXPECT_EQ(run_test(cfg).score(), 0);
    }

    {  // invalid band
        auto cfg_base = seqan3::align_cfg::mode{seqan3::global_alignment} |
                        seqan3::align_cfg::scoring{seqan3::nucleotide_scoring_scheme{}} |
                        seqan3::align_cfg::gap{seqan3::gap_scheme{seqan3::gap_score{-1}, seqan3::gap_open_score{-10}}};
        auto cfg_lower = cfg_base | seqan3::align_cfg::band{seqan3::static_band{seqan3::lower_bound{-10},
                                                                                seqan3::upper_bound{-5}}};
        auto cfg_upper = cfg_base | seqan3::align_cfg::band{seqan3::static_band{seqan3::lower_bound{5},
                                                                                seqan3::upper_bound{6}}};

        EXPECT_THROW(run_test(cfg_lower), seqan3::invalid_alignment_configuration);
        EXPECT_THROW(run_test(cfg_upper), seqan3::invalid_alignment_configuration);
    }
}

TEST(alignment_configurator, configure_affine_global_banded_with_alignment)
{
    auto cfg = seqan3::align_cfg::mode{seqan3::global_alignment} |
               seqan3::align_cfg::gap{seqan3::gap_scheme{seqan3::gap_score{-1}, seqan3::gap_open_score{-10}}} |
               seqan3::align_cfg::scoring{seqan3::nucleotide_scoring_scheme{}} |
               seqan3::align_cfg::band{seqan3::static_band{seqan3::lower_bound{-1}, seqan3::upper_bound{1}}};

    auto cfg_trace = cfg | seqan3::align_cfg::result{seqan3::with_alignment};
    auto cfg_begin = cfg | seqan3::align_cfg::result{seqan3::with_front_coordinate};
    auto cfg_end = cfg | seqan3::align_cfg::result{seqan3::with_back_coordinate};

    EXPECT_EQ(run_test(cfg_end).score(), 0);
    EXPECT_EQ(run_test(cfg_trace).score(), 0);
    EXPECT_EQ(run_test(cfg_begin).score(), 0);
}

TEST(alignment_configurator, configure_affine_global_semi)
{
    auto cfg = seqan3::align_cfg::mode{seqan3::global_alignment} |
               seqan3::align_cfg::scoring{seqan3::nucleotide_scoring_scheme{}} |
               seqan3::align_cfg::gap{seqan3::gap_scheme{seqan3::gap_score{-1}, seqan3::gap_open_score{-10}}} |
               seqan3::align_cfg::aligned_ends{seqan3::free_ends_all};

    EXPECT_EQ(run_test(cfg).score(), 0);
}

TEST(alignment_configurator, configure_affine_local)
{
    auto cfg = seqan3::align_cfg::mode{seqan3::local_alignment} |
               seqan3::align_cfg::gap{seqan3::gap_scheme{seqan3::gap_score{-1}, seqan3::gap_open_score{-10}}} |
               seqan3::align_cfg::scoring{seqan3::nucleotide_scoring_scheme{}};

    EXPECT_EQ(run_test(cfg).score(), 0);
}

TEST(alignment_configurator, configure_affine_local_back_coordinate)
{
    auto cfg = seqan3::align_cfg::mode{seqan3::local_alignment} |
               seqan3::align_cfg::gap{seqan3::gap_scheme{seqan3::gap_score{-1}, seqan3::gap_open_score{-10}}} |
               seqan3::align_cfg::scoring{seqan3::nucleotide_scoring_scheme{}} |
               seqan3::align_cfg::result{seqan3::with_back_coordinate};

    EXPECT_EQ(run_test(cfg).score(), 0);
}

TEST(alignment_configurator, configure_affine_local_front_coordinate)
{
    auto cfg = seqan3::align_cfg::mode{seqan3::local_alignment} |
               seqan3::align_cfg::gap{seqan3::gap_scheme{seqan3::gap_score{-1}, seqan3::gap_open_score{-10}}} |
               seqan3::align_cfg::scoring{seqan3::nucleotide_scoring_scheme{}} |
               seqan3::align_cfg::result{seqan3::with_front_coordinate};

    EXPECT_EQ(run_test(cfg).score(), 0);
}

TEST(alignment_configurator, configure_affine_local_alignment)
{
    auto cfg = seqan3::align_cfg::mode{seqan3::local_alignment} |
               seqan3::align_cfg::gap{seqan3::gap_scheme{seqan3::gap_score{-1}, seqan3::gap_open_score{-10}}} |
               seqan3::align_cfg::scoring{seqan3::nucleotide_scoring_scheme{}} |
               seqan3::align_cfg::result{seqan3::with_alignment};

    EXPECT_EQ(run_test(cfg).score(), 0);
}

TEST(alignment_configurator, configure_result_score_type)
{
    auto cfg = seqan3::align_cfg::edit | seqan3::align_cfg::result{seqan3::with_back_coordinate,
                                                                   seqan3::using_score_type<double>};
    auto result = run_test(cfg);

    EXPECT_DOUBLE_EQ(result.score(), 0.0);
    EXPECT_EQ(result.back_coordinate(),
              (seqan3::alignment_coordinate{seqan3::detail::column_index_type{4u},
                                            seqan3::detail::row_index_type{4u}}));
    EXPECT_TRUE((std::same_as<decltype(result.score()), double>));
}

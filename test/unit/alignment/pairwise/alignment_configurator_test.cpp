// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <range/v3/view/single.hpp>

#include <seqan3/alignment/pairwise/alignment_configurator.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/range/views/chunk.hpp>
#include <seqan3/range/views/zip.hpp>
#include <seqan3/std/ranges>

using namespace seqan3;

auto setup()
{
    auto data = std::tuple{"ACGT"_dna4, "ACGT"_dna4};
    return std::views::single(std::move(data));
}

template <typename config_t>
auto run_test(config_t const & cfg)
{
    auto r = setup();
    auto configuration_result = detail::alignment_configurator::configure<decltype(r)>(cfg);
    auto algorithm = configuration_result.first;

    auto indexed_sequence_pairs = views::zip(r, std::views::iota(0)) | views::chunk(1);
    return algorithm(*indexed_sequence_pairs.begin())[0];
}

TEST(alignment_configurator, configure_edit)
{
    EXPECT_EQ(run_test(align_cfg::edit).score(), 0);
}

TEST(alignment_configurator, configure_edit_end_position)
{
    EXPECT_EQ(run_test(align_cfg::edit | align_cfg::result{with_back_coordinate}).score(), 0);
}

TEST(alignment_configurator, configure_edit_begin_position)
{
    EXPECT_EQ(run_test(align_cfg::edit | align_cfg::result{with_front_coordinate}).score(), 0);
}

TEST(alignment_configurator, configure_edit_trace)
{
    EXPECT_EQ(run_test(align_cfg::edit | align_cfg::result{with_alignment}).score(), 0);
}

TEST(alignment_configurator, configure_edit_semi)
{
    EXPECT_EQ(run_test(align_cfg::edit | align_cfg::aligned_ends{free_ends_first}).score(), 0);
}

TEST(alignment_configurator, configure_edit_banded)
{
    EXPECT_THROW((run_test(align_cfg::edit | align_cfg::band{static_band{lower_bound{-1}, upper_bound{1}}})),
                 invalid_alignment_configuration);
}

TEST(alignment_configurator, configure_edit_max_error)
{
    EXPECT_EQ(run_test(align_cfg::edit | align_cfg::max_error{3u}).score(), 0);
}

TEST(alignment_configurator, configure_affine_global)
{
    auto cfg = align_cfg::mode{global_alignment} |
               align_cfg::gap{gap_scheme{gap_score{-1}, gap_open_score{-10}}} |
               align_cfg::scoring{nucleotide_scoring_scheme{}};

    EXPECT_EQ(run_test(cfg).score(), 0);
}

TEST(alignment_configurator, configure_affine_global_max_error)
{
    auto cfg = align_cfg::mode{global_alignment} |
               align_cfg::gap{gap_scheme{gap_score{-1}, gap_open_score{-10}}} |
               align_cfg::scoring{nucleotide_scoring_scheme{}} |
               align_cfg::max_error{5u};

    EXPECT_THROW(run_test(cfg), invalid_alignment_configuration);
}

TEST(alignment_configurator, configure_affine_global_end_position)
{
    auto cfg = align_cfg::mode{global_alignment} |
               align_cfg::gap{gap_scheme{gap_score{-1}, gap_open_score{-10}}} |
               align_cfg::scoring{nucleotide_scoring_scheme{}} |
               align_cfg::result{with_back_coordinate};

    EXPECT_EQ(run_test(cfg).score(), 0);
}

TEST(alignment_configurator, configure_affine_global_begin_position)
{
    auto cfg = align_cfg::mode{global_alignment} |
               align_cfg::gap{gap_scheme{gap_score{-1}, gap_open_score{-10}}} |
               align_cfg::scoring{nucleotide_scoring_scheme{}} |
               align_cfg::result{with_front_coordinate};

    EXPECT_EQ(run_test(cfg).score(), 0);
}

TEST(alignment_configurator, configure_affine_global_trace)
{
    auto cfg = align_cfg::mode{global_alignment} |
               align_cfg::gap{gap_scheme{gap_score{-1}, gap_open_score{-10}}} |
               align_cfg::scoring{nucleotide_scoring_scheme{}} |
               align_cfg::result{with_alignment};

    EXPECT_EQ(run_test(cfg).score(), 0);
}

TEST(alignment_configurator, configure_affine_global_banded)
{
    {
        auto cfg = align_cfg::mode{global_alignment} |
                   align_cfg::scoring{nucleotide_scoring_scheme{}} |
                   align_cfg::gap{gap_scheme{gap_score{-1}, gap_open_score{-10}}} |
                   align_cfg::band{static_band{lower_bound{-1}, upper_bound{1}}};

        EXPECT_EQ(run_test(cfg).score(), 0);
    }

    {  // invalid band
        auto cfg_base = align_cfg::mode{global_alignment} |
                        align_cfg::scoring{nucleotide_scoring_scheme{}} |
                        align_cfg::gap{gap_scheme{gap_score{-1}, gap_open_score{-10}}};
        auto cfg_lower = cfg_base | align_cfg::band{static_band{lower_bound{-10}, upper_bound{-5}}};
        auto cfg_upper = cfg_base | align_cfg::band{static_band{lower_bound{5}, upper_bound{6}}};

        EXPECT_THROW(run_test(cfg_lower), invalid_alignment_configuration);
        EXPECT_THROW(run_test(cfg_upper), invalid_alignment_configuration);
    }
}

TEST(alignment_configurator, configure_affine_global_banded_with_alignment)
{
    auto cfg = align_cfg::mode{global_alignment} |
               align_cfg::gap{gap_scheme{gap_score{-1}, gap_open_score{-10}}} |
               align_cfg::scoring{nucleotide_scoring_scheme{}} |
               align_cfg::band{static_band{lower_bound{-1}, upper_bound{1}}};

    auto cfg_trace = cfg | align_cfg::result{with_alignment};
    auto cfg_begin = cfg | align_cfg::result{with_front_coordinate};
    auto cfg_end = cfg | align_cfg::result{with_back_coordinate};

    EXPECT_EQ(run_test(cfg_end).score(), 0);
    EXPECT_EQ(run_test(cfg_trace).score(), 0);
    EXPECT_EQ(run_test(cfg_begin).score(), 0);
}

TEST(alignment_configurator, configure_affine_global_semi)
{
    auto cfg = align_cfg::mode{global_alignment} |
               align_cfg::scoring{nucleotide_scoring_scheme{}} |
               align_cfg::gap{gap_scheme{gap_score{-1}, gap_open_score{-10}}} |
               align_cfg::aligned_ends{free_ends_all};

    EXPECT_EQ(run_test(cfg).score(), 0);
}

TEST(alignment_configurator, configure_affine_local)
{
    auto cfg = align_cfg::mode{local_alignment} |
               align_cfg::gap{gap_scheme{gap_score{-1}, gap_open_score{-10}}} |
               align_cfg::scoring{nucleotide_scoring_scheme{}};

    EXPECT_EQ(run_test(cfg).score(), 0);
}

TEST(alignment_configurator, configure_affine_local_back_coordinate)
{
    auto cfg = align_cfg::mode{local_alignment} |
               align_cfg::gap{gap_scheme{gap_score{-1}, gap_open_score{-10}}} |
               align_cfg::scoring{nucleotide_scoring_scheme{}} |
               align_cfg::result{with_back_coordinate};

    EXPECT_EQ(run_test(cfg).score(), 0);
}

TEST(alignment_configurator, configure_affine_local_front_coordinate)
{
    auto cfg = align_cfg::mode{local_alignment} |
               align_cfg::gap{gap_scheme{gap_score{-1}, gap_open_score{-10}}} |
               align_cfg::scoring{nucleotide_scoring_scheme{}} |
               align_cfg::result{with_front_coordinate};

    EXPECT_EQ(run_test(cfg).score(), 0);
}

TEST(alignment_configurator, configure_affine_local_alignment)
{
    auto cfg = align_cfg::mode{local_alignment} |
               align_cfg::gap{gap_scheme{gap_score{-1}, gap_open_score{-10}}} |
               align_cfg::scoring{nucleotide_scoring_scheme{}} |
               align_cfg::result{with_alignment};

    EXPECT_EQ(run_test(cfg).score(), 0);
}

TEST(alignment_configurator, configure_result_score_type)
{
    auto cfg = align_cfg::edit | align_cfg::result{with_back_coordinate, using_score_type<double>};
    auto result = run_test(cfg);

    EXPECT_DOUBLE_EQ(result.score(), 0.0);
    EXPECT_EQ(result.back_coordinate(),
              (alignment_coordinate{detail::column_index_type{4u}, detail::row_index_type{4u}}));
    EXPECT_TRUE((std::same_as<decltype(result.score()), double>));
}

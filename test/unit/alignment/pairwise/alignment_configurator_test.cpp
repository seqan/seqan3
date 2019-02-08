// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <range/v3/view/single.hpp>

#include <seqan3/alignment/pairwise/alignment_configurator.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/range/view/persist.hpp>

using namespace seqan3;

auto setup()
{
    auto data = std::tuple{"ACGT"_dna4, "ACGT"_dna4};
    return ranges::view::single(std::move(data)) | view::persist;
}

TEST(alignment_configurator, configure_edit)
{
    auto cfg = align_cfg::edit;

    auto r = setup();
    auto fn = detail::alignment_configurator::configure(r, cfg);
    auto [seq1, seq2] = *seqan3::begin(r);

    EXPECT_EQ((fn(seq1, seq2).get_score()), 0);
}

TEST(alignment_configurator, configure_edit_end_position)
{
    auto cfg = align_cfg::edit | align_cfg::result{align_cfg::with_end_position};

    auto r = setup();
    auto fn = detail::alignment_configurator::configure(r, cfg);
    auto [seq1, seq2] = *seqan3::begin(r);

    EXPECT_EQ((fn(seq1, seq2).get_score()), 0);
}

TEST(alignment_configurator, configure_edit_begin_position)
{
    auto cfg = align_cfg::edit | align_cfg::result{align_cfg::with_begin_position};

    auto r = setup();
    auto fn = detail::alignment_configurator::configure(r, cfg);
    auto [seq1, seq2] = *seqan3::begin(r);

    EXPECT_EQ((fn(seq1, seq2).get_score()), 0);
}

TEST(alignment_configurator, configure_edit_trace)
{
    auto cfg = align_cfg::edit | align_cfg::result{align_cfg::with_trace};

    auto r = setup();
    auto fn = detail::alignment_configurator::configure(r, cfg);
    auto [seq1, seq2] = *seqan3::begin(r);

    EXPECT_EQ((fn(seq1, seq2).get_score()), 0);
}

TEST(alignment_configurator, configure_edit_semi)
{
    auto cfg = align_cfg::edit | align_cfg::aligned_ends{align_cfg::seq1_ends_free};
    EXPECT_THROW(detail::alignment_configurator::configure(setup(), cfg), invalid_alignment_configuration);
}

TEST(alignment_configurator, configure_edit_banded)
{
    auto cfg = align_cfg::edit | align_cfg::band{static_band{lower_bound{-1}, upper_bound{1}}};
    EXPECT_THROW(detail::alignment_configurator::configure(setup(), cfg), invalid_alignment_configuration);
}

TEST(alignment_configurator, configure_affine_global)
{
    auto cfg = align_cfg::mode{align_cfg::global_alignment} |
               align_cfg::gap{gap_scheme{gap_score{-1}, gap_open_score{-10}}} |
               align_cfg::scoring{nucleotide_scoring_scheme{}};

    auto r = setup();
    auto fn = detail::alignment_configurator::configure(r, cfg);
    auto [seq1, seq2] = *seqan3::begin(r);

    EXPECT_EQ((fn(seq1, seq2).get_score()), 0);
}

TEST(alignment_configurator, configure_affine_global_end_position)
{
    auto cfg = align_cfg::mode{align_cfg::global_alignment} |
               align_cfg::gap{gap_scheme{gap_score{-1}, gap_open_score{-10}}} |
               align_cfg::scoring{nucleotide_scoring_scheme{}} |
               align_cfg::result{align_cfg::with_end_position};

    auto r = setup();
    auto fn = detail::alignment_configurator::configure(r, cfg);
    auto [seq1, seq2] = *seqan3::begin(r);

    EXPECT_EQ((fn(seq1, seq2).get_score()), 0);
}

TEST(alignment_configurator, configure_affine_global_begin_position)
{
    auto cfg = align_cfg::mode{align_cfg::global_alignment} |
               align_cfg::gap{gap_scheme{gap_score{-1}, gap_open_score{-10}}} |
               align_cfg::scoring{nucleotide_scoring_scheme{}} |
               align_cfg::result{align_cfg::with_begin_position};

    EXPECT_THROW(detail::alignment_configurator::configure(setup(), cfg), invalid_alignment_configuration);
}

TEST(alignment_configurator, configure_affine_global_trace)
{
    auto cfg = align_cfg::mode{align_cfg::global_alignment} |
               align_cfg::gap{gap_scheme{gap_score{-1}, gap_open_score{-10}}} |
               align_cfg::scoring{nucleotide_scoring_scheme{}} |
               align_cfg::result{align_cfg::with_trace};

    EXPECT_THROW(detail::alignment_configurator::configure(setup(), cfg), invalid_alignment_configuration);
}

TEST(alignment_configurator, configure_affine_global_banded)
{
    {
        auto cfg = align_cfg::mode{align_cfg::global_alignment} |
                   align_cfg::scoring{nucleotide_scoring_scheme{}} |
                   align_cfg::gap{gap_scheme{gap_score{-1}, gap_open_score{-10}}} |
                   align_cfg::band{static_band{lower_bound{-1}, upper_bound{1}}};

        auto r = setup();
        auto fn = detail::alignment_configurator::configure(r, cfg);
        auto [seq1, seq2] = *seqan3::begin(r);

        EXPECT_EQ((fn(seq1, seq2).get_score()), 0);
    }

    {  // invalid band
        auto cfg_base = align_cfg::mode{align_cfg::global_alignment} |
                        align_cfg::scoring{nucleotide_scoring_scheme{}} |
                        align_cfg::gap{gap_scheme{gap_score{-1}, gap_open_score{-10}}};
        auto cfg_lower = cfg_base | align_cfg::band{static_band{lower_bound{-10}, upper_bound{-5}}};
        auto cfg_upper = cfg_base | align_cfg::band{static_band{lower_bound{5}, upper_bound{6}}};

        auto r = setup();

        {
            auto fn = detail::alignment_configurator::configure(r, cfg_lower);
            auto [seq1, seq2] = *seqan3::begin(r);

            EXPECT_THROW((fn(seq1, seq2).get_score()), invalid_alignment_configuration);
        }

        {
            auto fn = detail::alignment_configurator::configure(r, cfg_upper);
            auto [seq1, seq2] = *seqan3::begin(r);

            EXPECT_THROW((fn(seq1, seq2).get_score()), invalid_alignment_configuration);
        }
    }
}

TEST(alignment_configurator, configure_affine_global_semi)
{
    auto cfg = align_cfg::mode{align_cfg::global_alignment} |
               align_cfg::scoring{nucleotide_scoring_scheme{}} |
               align_cfg::gap{gap_scheme{gap_score{-1}, gap_open_score{-10}}} |
               align_cfg::aligned_ends{align_cfg::all_ends_free};

    auto r = setup();
    auto fn = detail::alignment_configurator::configure(r, cfg);
    auto [seq1, seq2] = *seqan3::begin(r);

    EXPECT_EQ((fn(seq1, seq2).get_score()), 0);
}

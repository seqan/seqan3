// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alignment/configuration/align_config_scoring.hpp>
#include <seqan3/alignment/scoring/aminoacid_scoring_scheme.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/core/algorithm/configuration.hpp>
#include <seqan3/std/concepts>

template <typename t>
class align_confg_scoring_test : public ::testing::Test
{
public:

    using scheme_t = std::tuple_element_t<0, t>;
    using alph_t   = std::tuple_element_t<1, t>;
};

using test_types = ::testing::Types<std::tuple<seqan3::aminoacid_scoring_scheme<int8_t>, seqan3::aa27>,
                                    std::tuple<seqan3::nucleotide_scoring_scheme<int8_t>, seqan3::dna15>>;
TYPED_TEST_SUITE(align_confg_scoring_test, test_types, );

TYPED_TEST(align_confg_scoring_test, config_element)
{
    using scheme_t = typename TestFixture::scheme_t;
    EXPECT_TRUE((seqan3::detail::config_element<seqan3::align_cfg::scoring<scheme_t>>));
}

TYPED_TEST(align_confg_scoring_test, configuration)
{
    using alph_t   = typename TestFixture::alph_t;
    using scheme_t = typename TestFixture::scheme_t;
    {
        seqan3::align_cfg::scoring elem{scheme_t{}};
        seqan3::configuration cfg{elem};

        EXPECT_EQ((seqan3::get<seqan3::align_cfg::scoring>(cfg).value.score(seqan3::assign_char_to('a', alph_t{}),
                                                                            seqan3::assign_char_to('a', alph_t{}))), 0);
    }

    {
        seqan3::configuration cfg{seqan3::align_cfg::scoring{scheme_t{}}};

        EXPECT_EQ((seqan3::get<seqan3::align_cfg::scoring>(cfg).value.score(seqan3::assign_char_to('a', alph_t{}),
                                                                            seqan3::assign_char_to('c', alph_t{}))), -1);
    }
}

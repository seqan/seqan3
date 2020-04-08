// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

#include "fm_index_collection_test_template.hpp"
#include "fm_index_test_template.hpp"

using t1 = std::pair<seqan3::fm_index<seqan3::dna4, seqan3::text_layout::single>, seqan3::dna4_vector>;
INSTANTIATE_TYPED_TEST_SUITE_P(dna4, fm_index_test, t1, );
using t2 = std::pair<seqan3::fm_index<seqan3::dna4, seqan3::text_layout::collection>, std::vector<seqan3::dna4_vector>>;
INSTANTIATE_TYPED_TEST_SUITE_P(dna4_collection, fm_index_collection_test, t2, );

TEST(fm_index_test, additional_concepts)
{
    EXPECT_TRUE(seqan3::detail::sdsl_index<seqan3::default_sdsl_index_type>);
}

TEST(fm_index_test, cerealisation_errors)
{
#if SEQAN3_WITH_CEREAL

    using seqan3::operator""_dna4;

    seqan3::fm_index<seqan3::dna4, seqan3::text_layout::single> index{"AGTCTGATGCTGCTAC"_dna4};

    seqan3::test::tmp_filename filename{"cereal_test"};

    {
        std::ofstream os{filename.get_path(), std::ios::binary};
        cereal::BinaryOutputArchive oarchive{os};
        oarchive(index);
    }

    {
        seqan3::fm_index<seqan3::dna5, seqan3::text_layout::single> in;
        std::ifstream is{filename.get_path(), std::ios::binary};
        cereal::BinaryInputArchive iarchive{is};
        EXPECT_THROW(iarchive(in), std::logic_error);
    }

    {
        seqan3::fm_index<seqan3::dna4, seqan3::text_layout::collection> in;
        std::ifstream is{filename.get_path(), std::ios::binary};
        cereal::BinaryInputArchive iarchive{is};
        EXPECT_THROW(iarchive(in), std::logic_error);
    }
#endif
}

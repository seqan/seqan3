// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include "fm_index_collection_test_template.hpp"
#include "fm_index_test_template.hpp"

using t1 = std::pair<fm_index<dna4, text_layout::single>, std::vector<dna4>>;
INSTANTIATE_TYPED_TEST_SUITE_P(dna4, fm_index_test, t1, );
using t2 = std::pair<fm_index<dna4, text_layout::collection>, std::vector<std::vector<dna4>>>;
INSTANTIATE_TYPED_TEST_SUITE_P(dna4_collection, fm_index_collection_test, t2, );

TEST(fm_index_test, additional_concepts)
{
    EXPECT_TRUE(detail::sdsl_index<default_sdsl_index_type>);
}

TEST(fm_index_test, cerealisation_errors)
{
#if SEQAN3_WITH_CEREAL
    using namespace seqan3;

    fm_index<dna4, text_layout::single> index{"AGTCTGATGCTGCTAC"_dna4};

    test::tmp_filename filename{"cereal_test"};

    {
        std::ofstream os{filename.get_path(), std::ios::binary};
        cereal::BinaryOutputArchive oarchive{os};
        oarchive(index);
    }

    {
        fm_index<dna5, text_layout::single> in;
        std::ifstream is{filename.get_path(), std::ios::binary};
        cereal::BinaryInputArchive iarchive{is};
        EXPECT_THROW(iarchive(in), std::logic_error);
    }

    {
        fm_index<dna4, text_layout::collection> in;
        std::ifstream is{filename.get_path(), std::ios::binary};
        cereal::BinaryInputArchive iarchive{is};
        EXPECT_THROW(iarchive(in), std::logic_error);
    }
#endif
}

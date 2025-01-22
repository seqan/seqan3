// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/test/tmp_directory.hpp>

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

    seqan3::test::tmp_directory tmp;
    auto filename = tmp.path() / "cereal_test";

    {
        std::ofstream os{filename, std::ios::binary};
        cereal::BinaryOutputArchive oarchive{os};
        oarchive(index);
    }

    {
        seqan3::fm_index<seqan3::dna5, seqan3::text_layout::single> in;
        std::ifstream is{filename, std::ios::binary};
        cereal::BinaryInputArchive iarchive{is};
        EXPECT_THROW(iarchive(in), std::logic_error);
    }

    {
        seqan3::fm_index<seqan3::dna4, seqan3::text_layout::collection> in;
        std::ifstream is{filename, std::ios::binary};
        cereal::BinaryInputArchive iarchive{is};
        EXPECT_THROW(iarchive(in), std::logic_error);
    }
#endif
}

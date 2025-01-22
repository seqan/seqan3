// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/alphabet/detail/debug_stream_alphabet.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/io/detail/record_like.hpp>
#include <seqan3/io/sequence_file/record.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/test/expect_same_type.hpp>
#include <seqan3/utility/tuple/concept.hpp>
#include <seqan3/utility/type_list/type_list.hpp>

using seqan3::operator""_dna4;

// ----------------------------------------------------------------------------
// record
// ----------------------------------------------------------------------------

struct sequence_record : public ::testing::Test
{
    using types = seqan3::type_list<std::string, seqan3::dna4_vector>;
    using types_as_ids = seqan3::fields<seqan3::field::id, seqan3::field::seq>;
    using record_type = seqan3::sequence_record<types, types_as_ids>;
};

TEST_F(sequence_record, concept)
{
    EXPECT_TRUE((seqan3::detail::record_like<record_type>));
}

TEST_F(sequence_record, definition_tuple_traits)
{
    EXPECT_SAME_TYPE((std::tuple<std::string, seqan3::dna4_vector>), typename record_type::base_type);

    EXPECT_SAME_TYPE(std::string, (std::tuple_element_t<0, record_type>));
    EXPECT_SAME_TYPE(seqan3::dna4_vector, (std::tuple_element_t<1, record_type>));
    EXPECT_EQ(std::tuple_size_v<record_type>, 2ul);

    EXPECT_TRUE(seqan3::tuple_like<record_type>);
}

TEST_F(sequence_record, construction)
{
    [[maybe_unused]] record_type r{"MY ID", "ACGT"_dna4};
}

TEST_F(sequence_record, get_by_index)
{
    record_type r{"MY ID", "ACGT"_dna4};

    EXPECT_EQ(std::get<0>(r), "MY ID");
    EXPECT_RANGE_EQ(std::get<1>(r), "ACGT"_dna4);
}

TEST_F(sequence_record, get_by_type)
{
    record_type r{"MY ID", "ACGT"_dna4};

    EXPECT_EQ(std::get<std::string>(r), "MY ID");
    EXPECT_RANGE_EQ(std::get<seqan3::dna4_vector>(r), "ACGT"_dna4);
}

TEST_F(sequence_record, get_by_member)
{
    record_type r{"MY ID", "ACGT"_dna4};

    EXPECT_EQ(r.id(), "MY ID");
    EXPECT_RANGE_EQ(r.sequence(), "ACGT"_dna4);
}

TEST_F(sequence_record, member_types)
{
    record_type r{"MY ID", "ACGT"_dna4};

    EXPECT_SAME_TYPE(std::string &, decltype(r.id()));
    EXPECT_SAME_TYPE(seqan3::dna4_vector &, decltype(r.sequence()));

    EXPECT_SAME_TYPE(std::string const &, decltype(std::as_const(r).id()));
    EXPECT_SAME_TYPE(seqan3::dna4_vector const &, decltype(std::as_const(r).sequence()));

    EXPECT_SAME_TYPE(std::string &&, decltype(std::move(r).id()));
    EXPECT_SAME_TYPE(seqan3::dna4_vector &&, decltype(std::move(r).sequence()));

    EXPECT_SAME_TYPE(std::string const &&, decltype(std::move(std::as_const(r)).id()));
    EXPECT_SAME_TYPE(seqan3::dna4_vector const &&, decltype(std::move(std::as_const(r)).sequence()));
}

// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <seqan3/std/algorithm>
#include <sstream>

#include <gtest/gtest.h>

#include <seqan3/alphabet/detail/debug_stream_alphabet.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/io/record.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/test/expect_same_type.hpp>
#include <seqan3/utility/tuple/concept.hpp>

using seqan3::operator""_dna4;

using default_fields = seqan3::fields<seqan3::field::seq, seqan3::field::id, seqan3::field::qual>;

// This is needed for EXPECT_RANGE_EQ:
namespace seqan3
{
template <typename char_t>
inline debug_stream_type<char_t> & operator<<(debug_stream_type<char_t> & stream, field f)
{
    stream << "<field: " << static_cast<size_t>(f) << ">";
    return stream;
}
} // namespace seqan3

// ----------------------------------------------------------------------------
// fields
// ----------------------------------------------------------------------------

TEST(fields, usage)
{
    std::array comp{seqan3::field::seq, seqan3::field::id, seqan3::field::qual};

    EXPECT_RANGE_EQ(default_fields::as_array, comp);
    EXPECT_TRUE(default_fields::contains(seqan3::field::seq));
    EXPECT_TRUE(default_fields::contains(seqan3::field::id));
    EXPECT_TRUE(default_fields::contains(seqan3::field::qual));
#ifdef SEQAN3_DEPRECATED_310
    EXPECT_FALSE(default_fields::contains(seqan3::field::_seq_qual_deprecated));
#endif // SEQAN3_DEPRECATED_310
    EXPECT_EQ(default_fields::index_of(seqan3::field::seq), 0ul);
    EXPECT_EQ(default_fields::index_of(seqan3::field::id),  1ul);
    EXPECT_EQ(default_fields::index_of(seqan3::field::qual), 2ul);
}

// ----------------------------------------------------------------------------
// record
// ----------------------------------------------------------------------------

struct record : public ::testing::Test
{
    using types        = seqan3::type_list<std::string, seqan3::dna4_vector>;
    using types_as_ids = seqan3::fields<seqan3::field::id, seqan3::field::seq>;
    using record_type  = seqan3::record<types, types_as_ids>;
};

TEST_F(record, definition_tuple_traits)
{
    EXPECT_TRUE((std::is_same_v<typename record_type::base_type,
                                std::tuple<std::string, seqan3::dna4_vector>>));

    EXPECT_TRUE((std::is_same_v<std::tuple_element_t<0, record_type>,
                                std::string>));
    EXPECT_TRUE((std::is_same_v<std::tuple_element_t<1, record_type>,
                                seqan3::dna4_vector>));
    EXPECT_EQ(std::tuple_size_v<record_type>, 2ul);

    EXPECT_TRUE(seqan3::tuple_like<record_type>);
}

TEST_F(record, construction)
{
    [[maybe_unused]] record_type r{"MY ID", "ACGT"_dna4};
}

TEST_F(record, get_by_index)
{
    record_type r{"MY ID", "ACGT"_dna4};

    EXPECT_EQ(std::get<0>(r), "MY ID");
    EXPECT_RANGE_EQ(std::get<1>(r), "ACGT"_dna4);
}

TEST_F(record, get_by_type)
{
    record_type r{"MY ID", "ACGT"_dna4};

    EXPECT_EQ(std::get<std::string>(r), "MY ID");
    EXPECT_RANGE_EQ(std::get<seqan3::dna4_vector>(r), "ACGT"_dna4);
}

#ifdef SEQAN3_DEPRECATED_310
TEST_F(record, get_by_field)
{
    record_type r{"MY ID", "ACGT"_dna4};

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
    EXPECT_EQ(seqan3::get<seqan3::field::id>(r), "MY ID");
    EXPECT_RANGE_EQ(seqan3::get<seqan3::field::seq>(r), "ACGT"_dna4);
#pragma GCC diagnostic pop
}

TEST_F(record, gcc_issue_94967)
{
    record_type r{"MY ID", "ACGT"_dna4};

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
    // gcc 7 implemented a different std::get(const &&) definition which was later changed by
    // https://wg21.link/lwg2485
    EXPECT_SAME_TYPE(std::string const &&, decltype(seqan3::get<seqan3::field::id>(std::move(std::as_const(r)))));
#pragma GCC diagnostic pop
}
#endif // SEQAN3_DEPRECATED_310

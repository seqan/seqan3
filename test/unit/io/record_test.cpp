// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <sstream>

#include <gtest/gtest.h>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/io/record.hpp>
#include <seqan3/io/detail/record.hpp>
#include <seqan3/core/concept/tuple.hpp>
#include <seqan3/core/detail/reflection.hpp>
#include <seqan3/range/view/to_char.hpp>
#include <seqan3/std/algorithm>

using namespace seqan3;

// ----------------------------------------------------------------------------
// fields
// ----------------------------------------------------------------------------

TEST(fields, usage)
{
    using fields_t = fields<field::SEQ, field::ID, field::QUAL>;
    std::array comp{field::SEQ, field::ID, field::QUAL};

    EXPECT_TRUE(std::ranges::equal(fields_t::as_array, comp));
    EXPECT_TRUE(fields_t::contains(field::SEQ));
    EXPECT_TRUE(fields_t::contains(field::ID));
    EXPECT_TRUE(fields_t::contains(field::QUAL));
    EXPECT_FALSE(fields_t::contains(field::SEQ_QUAL));
    EXPECT_EQ(fields_t::index_of(field::SEQ), 0ul);
    EXPECT_EQ(fields_t::index_of(field::ID),  1ul);
    EXPECT_EQ(fields_t::index_of(field::QUAL), 2ul);
}

// ----------------------------------------------------------------------------
// record
// ----------------------------------------------------------------------------

struct record_ : public ::testing::Test
{
    using types        = type_list<std::string, dna4_vector>;
    using types_as_ids = fields<field::ID, field::SEQ>;
    using record_type  = record<types, types_as_ids>;
};

TEST_F(record_, definition_tuple_traits)
{
    EXPECT_TRUE((std::is_same_v<typename record_type::base_type,
                                std::tuple<std::string, dna4_vector>>));

    EXPECT_TRUE((std::is_same_v<std::tuple_element_t<0, record_type>,
                                std::string>));
    EXPECT_TRUE((std::is_same_v<std::tuple_element_t<1, record_type>,
                                dna4_vector>));
    EXPECT_EQ(std::tuple_size_v<record_type>, 2ul);

    EXPECT_TRUE(TupleLike<record_type>);
}

TEST_F(record_, construction)
{
    [[maybe_unused]] record_type r{"MY ID", "ACGT"_dna4};
}

TEST_F(record_, get_by_index)
{
    record_type r{"MY ID", "ACGT"_dna4};

    EXPECT_EQ(get<0>(r), "MY ID");
    EXPECT_TRUE(std::ranges::equal(get<1>(r), "ACGT"_dna4));
}

TEST_F(record_, get_by_type)
{
    record_type r{"MY ID", "ACGT"_dna4};

    EXPECT_EQ(get<std::string>(r), "MY ID");
    EXPECT_TRUE(std::ranges::equal(get<dna4_vector>(r), "ACGT"_dna4));
}

TEST_F(record_, get_by_field)
{
    record_type r{"MY ID", "ACGT"_dna4};

    EXPECT_EQ(get<field::ID>(r), "MY ID");
    EXPECT_TRUE(std::ranges::equal(get<field::SEQ>(r), "ACGT"_dna4));
}

// ----------------------------------------------------------------------------
// detail/record.hpp
// ----------------------------------------------------------------------------

TEST(detail, select_types_with_ids)
{
    using types         = type_list<std::string, dna4_vector, std::vector<phred42>>;
    using types_as_ids  = fields<field::ID, field::SEQ, field::QUAL>;
    using selected_ids  = fields<field::QUAL, field::ID>;

    using selected_types = typename detail::select_types_with_ids<types, types_as_ids, selected_ids>::type;

    EXPECT_TRUE((std::is_same_v<selected_types,
                                type_list<std::vector<phred42>, std::string>>));
}

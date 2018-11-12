// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
//
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================

#include <sstream>

#include <gtest/gtest.h>

#include <range/v3/algorithm/equal.hpp>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/io/record.hpp>
#include <seqan3/io/detail/record.hpp>
#include <seqan3/core/concept/tuple.hpp>
#include <seqan3/core/detail/reflection.hpp>
#include <seqan3/range/view/to_char.hpp>

using namespace seqan3;
using namespace seqan3::literal;

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

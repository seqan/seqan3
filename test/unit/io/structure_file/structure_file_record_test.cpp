// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/alphabet/detail/debug_stream_alphabet.hpp>
#include <seqan3/alphabet/nucleotide/rna5.hpp>
#include <seqan3/alphabet/structure/wuss.hpp>
#include <seqan3/io/detail/record_like.hpp>
#include <seqan3/io/structure_file/record.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/test/expect_same_type.hpp>
#include <seqan3/utility/tuple/concept.hpp>
#include <seqan3/utility/type_list/type_list.hpp>

using seqan3::operator""_rna5;
using seqan3::operator""_wuss51;

// ----------------------------------------------------------------------------
// record
// ----------------------------------------------------------------------------

struct structure_record : public ::testing::Test
{
    using types = seqan3::type_list<std::string, seqan3::rna5_vector, std::vector<seqan3::wuss51>, double>;
    using types_as_ids =
        seqan3::fields<seqan3::field::id, seqan3::field::seq, seqan3::field::structure, seqan3::field::energy>;
    using record_type = seqan3::structure_record<types, types_as_ids>;
};

TEST_F(structure_record, concept)
{
    EXPECT_TRUE((seqan3::detail::record_like<record_type>));
}

TEST_F(structure_record, definition_tuple_traits)
{
    EXPECT_SAME_TYPE((std::tuple<std::string, seqan3::rna5_vector, std::vector<seqan3::wuss51>, double>),
                     typename record_type::base_type);

    EXPECT_SAME_TYPE(std::string, (std::tuple_element_t<0, record_type>));
    EXPECT_SAME_TYPE(seqan3::rna5_vector, (std::tuple_element_t<1, record_type>));
    EXPECT_SAME_TYPE(std::vector<seqan3::wuss51>, (std::tuple_element_t<2, record_type>));
    EXPECT_SAME_TYPE(double, (std::tuple_element_t<3, record_type>));
    EXPECT_EQ(std::tuple_size_v<record_type>, 4ul);

    EXPECT_TRUE(seqan3::tuple_like<record_type>);
}

TEST_F(structure_record, construction)
{
    [[maybe_unused]] record_type r{"MY ID", "ACGU"_rna5, "(())"_wuss51, 1.5};
}

TEST_F(structure_record, get_by_index)
{
    record_type r{"MY ID", "ACGU"_rna5, "(())"_wuss51, 1.5};

    EXPECT_EQ(std::get<0>(r), "MY ID");
    EXPECT_RANGE_EQ(std::get<1>(r), "ACGU"_rna5);
    EXPECT_RANGE_EQ(std::get<2>(r), "(())"_wuss51);
    EXPECT_DOUBLE_EQ(std::get<3>(r), 1.5);
}

TEST_F(structure_record, get_by_type)
{
    record_type r{"MY ID", "ACGU"_rna5, "(())"_wuss51, 1.5};

    EXPECT_EQ(std::get<std::string>(r), "MY ID");
    EXPECT_RANGE_EQ(std::get<seqan3::rna5_vector>(r), "ACGU"_rna5);
    EXPECT_RANGE_EQ(std::get<std::vector<seqan3::wuss51>>(r), "(())"_wuss51);
    EXPECT_DOUBLE_EQ(std::get<double>(r), 1.5);
}

TEST_F(structure_record, get_by_member)
{
    record_type r{"MY ID", "ACGU"_rna5, "(())"_wuss51, 1.5};

    EXPECT_EQ(r.id(), "MY ID");
    EXPECT_RANGE_EQ(r.sequence(), "ACGU"_rna5);
    EXPECT_RANGE_EQ(r.sequence_structure(), "(())"_wuss51);
    EXPECT_DOUBLE_EQ(r.energy(), 1.5);
}

TEST_F(structure_record, member_types)
{
    record_type r{"MY ID", "ACGU"_rna5, "(())"_wuss51, 1.5};

    EXPECT_SAME_TYPE(std::string &, decltype(r.id()));
    EXPECT_SAME_TYPE(seqan3::rna5_vector &, decltype(r.sequence()));
    EXPECT_SAME_TYPE(std::vector<seqan3::wuss51> &, decltype(r.sequence_structure()));
    EXPECT_SAME_TYPE(double &, decltype(r.energy()));

    EXPECT_SAME_TYPE(std::string const &, decltype(std::as_const(r).id()));
    EXPECT_SAME_TYPE(seqan3::rna5_vector const &, decltype(std::as_const(r).sequence()));
    EXPECT_SAME_TYPE(std::vector<seqan3::wuss51> const &, decltype(std::as_const(r).sequence_structure()));
    EXPECT_SAME_TYPE(double const &, decltype(std::as_const(r).energy()));

    EXPECT_SAME_TYPE(std::string &&, decltype(std::move(r).id()));
    EXPECT_SAME_TYPE(seqan3::rna5_vector &&, decltype(std::move(r).sequence()));
    EXPECT_SAME_TYPE(std::vector<seqan3::wuss51> &&, decltype(std::move(r).sequence_structure()));
    EXPECT_SAME_TYPE(double &&, decltype(std::move(r).energy()));

    EXPECT_SAME_TYPE(std::string const &&, decltype(std::move(std::as_const(r)).id()));
    EXPECT_SAME_TYPE(seqan3::rna5_vector const &&, decltype(std::move(std::as_const(r)).sequence()));
    EXPECT_SAME_TYPE(std::vector<seqan3::wuss51> const &&, decltype(std::move(std::as_const(r)).sequence_structure()));
    EXPECT_SAME_TYPE(double const &&, decltype(std::move(std::as_const(r)).energy()));
}

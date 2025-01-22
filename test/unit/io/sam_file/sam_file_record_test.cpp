// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/alphabet/cigar/cigar.hpp>
#include <seqan3/alphabet/detail/debug_stream_alphabet.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/core/debug_stream/byte.hpp>
#include <seqan3/core/debug_stream/optional.hpp>
#include <seqan3/core/debug_stream/tuple.hpp>
#include <seqan3/core/debug_stream/variant.hpp>
#include <seqan3/io/detail/record_like.hpp>
#include <seqan3/io/sam_file/header.hpp>
#include <seqan3/io/sam_file/record.hpp>
#include <seqan3/io/sam_file/sam_flag.hpp>
#include <seqan3/io/sam_file/sam_tag_dictionary.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/test/expect_same_type.hpp>
#include <seqan3/test/literal/cigar_literal.hpp>
#include <seqan3/utility/tuple/concept.hpp>
#include <seqan3/utility/type_list/type_list.hpp>

using seqan3::operator""_dna5;
using seqan3::operator""_phred42;
using seqan3::test::operator""_cigar;

// ----------------------------------------------------------------------------
// record
// ----------------------------------------------------------------------------

struct sam_record : public ::testing::Test
{
    using mate_t = std::tuple<std::optional<int32_t>, std::optional<int32_t>, int32_t>;

    using types = seqan3::type_list<std::string,                  // seqan3::field::id,
                                    seqan3::dna5_vector,          // seqan3::field::seq,
                                    std::vector<seqan3::phred42>, // seqan3::field::qual,
                                    std::optional<int32_t>,       // seqan3::field::ref_id,
                                    std::optional<int32_t>,       // seqan3::field::ref_offset,
                                    // seqan3::field::header_ptr,
                                    seqan3::sam_file_header<std::vector<std::string>> *,
                                    seqan3::sam_flag,            // seqan3::field::flag,
                                    mate_t,                      // seqan3::field::mate,
                                    uint8_t,                     // seqan3::field::mapq,
                                    std::vector<seqan3::cigar>,  // seqan3::field::cigar,
                                    seqan3::sam_tag_dictionary>; // seqan3::field::tags

    using types_as_ids = seqan3::fields<seqan3::field::id,
                                        seqan3::field::seq,
                                        seqan3::field::qual,
                                        seqan3::field::ref_id,
                                        seqan3::field::ref_offset,
                                        seqan3::field::header_ptr,
                                        seqan3::field::flag,
                                        seqan3::field::mate,
                                        seqan3::field::mapq,
                                        seqan3::field::cigar,
                                        seqan3::field::tags>;

    using record_type = seqan3::sam_record<types, types_as_ids>;

    record_type construct()
    {
        return {/*.id =*/"MY ID",
                /*.sequence =*/"ACGT"_dna5,
                /*.base_qualities =*/"!##$"_phred42,
                /*.reference_id =*/0, // "ref"
                /*.reference_position =*/0,
                /*.header_ptr =*/nullptr,
                /*.flag =*/seqan3::sam_flag{41u},
                {/*.mate_reference_id =*/0, /*.mate_position =*/9, /*.template_length =*/300},
                /*.mapping_quality =*/61u,
                /*.cigar_sequence =*/"1S1M1D1M1I"_cigar,
                /*.tags =*/seqan3::sam_tag_dictionary{}};
    }
};

TEST_F(sam_record, concept)
{
    EXPECT_TRUE((seqan3::detail::record_like<record_type>));
}

TEST_F(sam_record, definition_tuple_traits)
{
    EXPECT_SAME_TYPE((std::tuple<std::string,
                                 seqan3::dna5_vector,
                                 std::vector<seqan3::phred42>,
                                 std::optional<int32_t>,
                                 std::optional<int32_t>,
                                 seqan3::sam_file_header<std::vector<std::string>> *,
                                 seqan3::sam_flag,
                                 mate_t,
                                 uint8_t,
                                 std::vector<seqan3::cigar>,
                                 seqan3::sam_tag_dictionary>),
                     typename record_type::base_type);

    EXPECT_SAME_TYPE(std::string, (std::tuple_element_t<0, record_type>));
    EXPECT_SAME_TYPE(seqan3::dna5_vector, (std::tuple_element_t<1, record_type>));
    EXPECT_SAME_TYPE(std::vector<seqan3::phred42>, (std::tuple_element_t<2, record_type>));
    EXPECT_SAME_TYPE(std::optional<int32_t>, (std::tuple_element_t<3, record_type>));
    EXPECT_SAME_TYPE(std::optional<int32_t>, (std::tuple_element_t<4, record_type>));
    EXPECT_SAME_TYPE(seqan3::sam_file_header<std::vector<std::string>> *, (std::tuple_element_t<5, record_type>));
    EXPECT_SAME_TYPE(seqan3::sam_flag, (std::tuple_element_t<6, record_type>));
    EXPECT_SAME_TYPE(mate_t, (std::tuple_element_t<7, record_type>));
    EXPECT_SAME_TYPE(uint8_t, (std::tuple_element_t<8, record_type>));
    EXPECT_SAME_TYPE(std::vector<seqan3::cigar>, (std::tuple_element_t<9, record_type>));
    EXPECT_SAME_TYPE(seqan3::sam_tag_dictionary, (std::tuple_element_t<10, record_type>));

    EXPECT_EQ(std::tuple_size_v<record_type>, 11ul);

    EXPECT_TRUE(seqan3::tuple_like<record_type>);
}

TEST_F(sam_record, construction)
{
    [[maybe_unused]] record_type r{construct()};
}

TEST_F(sam_record, get_by_index)
{
    record_type r{construct()};

    EXPECT_EQ(std::get<0>(r), "MY ID");
    EXPECT_RANGE_EQ(std::get<1>(r), "ACGT"_dna5);
    EXPECT_RANGE_EQ(std::get<2>(r), "!##$"_phred42);
    EXPECT_EQ(std::get<3>(r), 0); // "ref"
    EXPECT_EQ(std::get<4>(r), 0);
    EXPECT_EQ(std::get<5>(r), nullptr);
    EXPECT_EQ(std::get<6>(r), seqan3::sam_flag{41u});
    EXPECT_EQ(std::get<7>(r), (mate_t{0, 9, 300}));
    EXPECT_EQ(std::get<8>(r), 61u);
    EXPECT_RANGE_EQ(std::get<9>(r), "1S1M1D1M1I"_cigar);
    EXPECT_EQ(std::get<10>(r), seqan3::sam_tag_dictionary{});
}

TEST_F(sam_record, get_by_type)
{
    record_type r{construct()};

    EXPECT_EQ(std::get<std::string>(r), "MY ID");
    EXPECT_RANGE_EQ(std::get<seqan3::dna5_vector>(r), "ACGT"_dna5);
    EXPECT_RANGE_EQ(std::get<std::vector<seqan3::phred42>>(r), "!##$"_phred42);
    // EXPECT_EQ(std::get<std::optional<int32_t>>(r), 0); // "ref" // ambiguous
    // EXPECT_EQ(std::get<std::optional<int32_t>>(r), 0); // ambiguous
    EXPECT_EQ(std::get<seqan3::sam_file_header<std::vector<std::string>> *>(r), nullptr);
    EXPECT_EQ(std::get<seqan3::sam_flag>(r), seqan3::sam_flag{41u});
    EXPECT_EQ(std::get<mate_t>(r), (mate_t{0, 9, 300}));
    EXPECT_EQ(std::get<uint8_t>(r), 61u);
    EXPECT_RANGE_EQ(std::get<std::vector<seqan3::cigar>>(r), "1S1M1D1M1I"_cigar);
    EXPECT_EQ(std::get<seqan3::sam_tag_dictionary>(r), seqan3::sam_tag_dictionary{});
}

TEST_F(sam_record, get_by_member)
{
    record_type r{construct()};

    EXPECT_EQ(r.id(), "MY ID");
    EXPECT_RANGE_EQ(r.sequence(), "ACGT"_dna5);
    EXPECT_RANGE_EQ(r.base_qualities(), "!##$"_phred42);
    EXPECT_EQ(r.reference_id(), 0); // "ref"
    EXPECT_EQ(r.reference_position(), 0);
    EXPECT_EQ(r.header_ptr(), nullptr);
    EXPECT_EQ(r.flag(), seqan3::sam_flag{41u});
    EXPECT_EQ(r.mate_reference_id(), 0);
    EXPECT_EQ(r.mate_position(), 9);
    EXPECT_EQ(r.template_length(), 300);
    EXPECT_EQ(r.mapping_quality(), 61u);
    EXPECT_RANGE_EQ(r.cigar_sequence(), "1S1M1D1M1I"_cigar);
    EXPECT_EQ(r.tags(), seqan3::sam_tag_dictionary{});
}

TEST_F(sam_record, member_types)
{
    record_type r{construct()};

    EXPECT_SAME_TYPE(std::string &, decltype(r.id()));
    EXPECT_SAME_TYPE(seqan3::dna5_vector &, decltype(r.sequence()));
    EXPECT_SAME_TYPE(std::vector<seqan3::phred42> &, decltype(r.base_qualities()));
    EXPECT_SAME_TYPE(std::optional<int32_t> &, decltype(r.reference_id()));
    EXPECT_SAME_TYPE(std::optional<int32_t> &, decltype(r.reference_position()));
    EXPECT_SAME_TYPE(seqan3::sam_file_header<std::vector<std::string>> *&, decltype(r.header_ptr()));
    EXPECT_SAME_TYPE(seqan3::sam_flag &, decltype(r.flag()));
    EXPECT_SAME_TYPE(std::optional<int32_t> &, decltype(r.mate_reference_id()));
    EXPECT_SAME_TYPE(std::optional<int32_t> &, decltype(r.mate_position()));
    EXPECT_SAME_TYPE(int32_t &, decltype(r.template_length()));
    EXPECT_SAME_TYPE(uint8_t &, decltype(r.mapping_quality()));
    EXPECT_SAME_TYPE(std::vector<seqan3::cigar> &, decltype(r.cigar_sequence()));
    EXPECT_SAME_TYPE(seqan3::sam_tag_dictionary &, decltype(r.tags()));

    EXPECT_SAME_TYPE(std::string const &, decltype(std::as_const(r.id())));
    EXPECT_SAME_TYPE(seqan3::dna5_vector const &, decltype(std::as_const(r.sequence())));
    EXPECT_SAME_TYPE(std::vector<seqan3::phred42> const &, decltype(std::as_const(r.base_qualities())));
    EXPECT_SAME_TYPE(std::optional<int32_t> const &, decltype(std::as_const(r.reference_id())));
    EXPECT_SAME_TYPE(std::optional<int32_t> const &, decltype(std::as_const(r.reference_position())));
    EXPECT_SAME_TYPE(seqan3::sam_file_header<std::vector<std::string>> * const &,
                     decltype(std::as_const(r.header_ptr())));
    EXPECT_SAME_TYPE(seqan3::sam_flag const &, decltype(std::as_const(r.flag())));
    EXPECT_SAME_TYPE(std::optional<int32_t> const &, decltype(std::as_const(r.mate_reference_id())));
    EXPECT_SAME_TYPE(std::optional<int32_t> const &, decltype(std::as_const(r.mate_position())));
    EXPECT_SAME_TYPE(int32_t const &, decltype(std::as_const(r.template_length())));
    EXPECT_SAME_TYPE(uint8_t const &, decltype(std::as_const(r.mapping_quality())));
    EXPECT_SAME_TYPE(std::vector<seqan3::cigar> const &, decltype(std::as_const(r.cigar_sequence())));
    EXPECT_SAME_TYPE(seqan3::sam_tag_dictionary const &, decltype(std::as_const(r.tags())));

    EXPECT_SAME_TYPE(std::string &&, decltype(std::move(r.id())));
    EXPECT_SAME_TYPE(seqan3::dna5_vector &&, decltype(std::move(r.sequence())));
    EXPECT_SAME_TYPE(std::vector<seqan3::phred42> &&, decltype(std::move(r.base_qualities())));
    EXPECT_SAME_TYPE(std::optional<int32_t> &&, decltype(std::move(r.reference_id())));
    EXPECT_SAME_TYPE(std::optional<int32_t> &&, decltype(std::move(r.reference_position())));
    EXPECT_SAME_TYPE(seqan3::sam_file_header<std::vector<std::string>> * &&, decltype(std::move(r.header_ptr())));
    EXPECT_SAME_TYPE(seqan3::sam_flag &&, decltype(std::move(r.flag())));
    EXPECT_SAME_TYPE(std::optional<int32_t> &&, decltype(std::move(r.mate_reference_id())));
    EXPECT_SAME_TYPE(std::optional<int32_t> &&, decltype(std::move(r.mate_position())));
    EXPECT_SAME_TYPE(int32_t &&, decltype(std::move(r.template_length())));
    EXPECT_SAME_TYPE(uint8_t &&, decltype(std::move(r.mapping_quality())));
    EXPECT_SAME_TYPE(std::vector<seqan3::cigar> &&, decltype(std::move(r.cigar_sequence())));
    EXPECT_SAME_TYPE(seqan3::sam_tag_dictionary &&, decltype(std::move(r.tags())));

    EXPECT_SAME_TYPE(std::string const &&, decltype(std::move(std::as_const(r.id()))));
    EXPECT_SAME_TYPE(seqan3::dna5_vector const &&, decltype(std::move(std::as_const(r.sequence()))));
    EXPECT_SAME_TYPE(std::vector<seqan3::phred42> const &&, decltype(std::move(std::as_const(r.base_qualities()))));
    EXPECT_SAME_TYPE(std::optional<int32_t> const &&, decltype(std::move(std::as_const(r.reference_id()))));
    EXPECT_SAME_TYPE(std::optional<int32_t> const &&, decltype(std::move(std::as_const(r.reference_position()))));
    EXPECT_SAME_TYPE(seqan3::sam_file_header<std::vector<std::string>> * const &&,
                     decltype(std::move(std::as_const(r.header_ptr()))));
    EXPECT_SAME_TYPE(seqan3::sam_flag const &&, decltype(std::move(std::as_const(r.flag()))));
    EXPECT_SAME_TYPE(std::optional<int32_t> const &&, decltype(std::move(std::as_const(r.mate_reference_id()))));
    EXPECT_SAME_TYPE(std::optional<int32_t> const &&, decltype(std::move(std::as_const(r.mate_position()))));
    EXPECT_SAME_TYPE(int32_t const &&, decltype(std::move(std::as_const(r.template_length()))));
    EXPECT_SAME_TYPE(uint8_t const &&, decltype(std::move(std::as_const(r.mapping_quality()))));
    EXPECT_SAME_TYPE(std::vector<seqan3::cigar> const &&, decltype(std::move(std::as_const(r.cigar_sequence()))));
    EXPECT_SAME_TYPE(seqan3::sam_tag_dictionary const &&, decltype(std::move(std::as_const(r.tags()))));
}

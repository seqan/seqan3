// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alphabet/cigar/cigar.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/core/detail/debug_stream_alphabet.hpp>
#include <seqan3/core/detail/debug_stream_byte.hpp>
#include <seqan3/core/detail/debug_stream_optional.hpp>
#include <seqan3/core/detail/debug_stream_tuple.hpp>
#include <seqan3/core/detail/debug_stream_variant.hpp>
#include <seqan3/io/detail/record_like.hpp>
#include <seqan3/io/sam_file/header.hpp>
#include <seqan3/io/sam_file/record.hpp>
#include <seqan3/io/sam_file/sam_flag.hpp>
#include <seqan3/io/sam_file/sam_tag_dictionary.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/test/expect_same_type.hpp>
#include <seqan3/utility/tuple/concept.hpp>
#include <seqan3/utility/type_list/type_list.hpp>

using seqan3::operator""_dna5;
using seqan3::operator""_phred42;
using seqan3::operator""_cigar_operation;

// ----------------------------------------------------------------------------
// record
// ----------------------------------------------------------------------------

struct sam_record : public ::testing::Test
{
    using alignment_t = std::pair<std::vector<seqan3::gapped<seqan3::dna5>>,
                                  std::vector<seqan3::gapped<seqan3::dna5>>>;

    using mate_t = std::tuple<std::optional<int32_t>, std::optional<int32_t>, int32_t>;

    using types = seqan3::type_list<std::string, // seqan3::field::id,
                                    seqan3::dna5_vector, // seqan3::field::seq,
                                    std::vector<seqan3::phred42>, // seqan3::field::qual,
                                    int32_t, // seqan3::field::offset,
                                    alignment_t, // seqan3::field::alignment,
                                    std::string, // seqan3::field::ref_id,
                                    std::optional<int32_t>, // seqan3::field::ref_offset,
                                    // seqan3::field::header_ptr,
                                    seqan3::sam_file_header<std::vector<std::string>> *,
                                    seqan3::sam_flag, // seqan3::field::flag,
                                    mate_t, // seqan3::field::mate,
                                    uint8_t, // seqan3::field::mapq,
                                    std::vector<seqan3::cigar>, // seqan3::field::cigar,
                                    seqan3::sam_tag_dictionary>; // seqan3::field::tags

    using types_as_ids = seqan3::fields<seqan3::field::id,
                                        seqan3::field::seq,
                                        seqan3::field::qual,
                                        seqan3::field::offset,
                                        seqan3::field::alignment,
                                        seqan3::field::ref_id,
                                        seqan3::field::ref_offset,
                                        seqan3::field::header_ptr,
                                        seqan3::field::flag,
                                        seqan3::field::mate,
                                        seqan3::field::mapq,
                                        seqan3::field::cigar,
                                        seqan3::field::tags>;

    using record_type  = seqan3::sam_record<types, types_as_ids>;

    record_type construct()
    {
        return
        {
            "MY ID",
            "ACGT"_dna5,
            "!##$"_phred42,
            1,
            alignment_t{},
            "ref",
            0,
            nullptr,
            seqan3::sam_flag{41u},
            {0, 9, 300},
            61u,
            {{1, 'S'_cigar_operation}, {1, 'M'_cigar_operation}, {1, 'D'_cigar_operation}, {1, 'M'_cigar_operation},
             {1, 'I'_cigar_operation}},
            seqan3::sam_tag_dictionary{}
        };
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
                                 int32_t,
                                 alignment_t,
                                 std::string,
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
    EXPECT_SAME_TYPE(int32_t, (std::tuple_element_t<3, record_type>));
    EXPECT_SAME_TYPE(alignment_t, (std::tuple_element_t<4, record_type>));
    EXPECT_SAME_TYPE(std::string, (std::tuple_element_t<5, record_type>));
    EXPECT_SAME_TYPE(std::optional<int32_t>, (std::tuple_element_t<6, record_type>));
    EXPECT_SAME_TYPE(seqan3::sam_file_header<std::vector<std::string>> *, (std::tuple_element_t<7, record_type>));
    EXPECT_SAME_TYPE(seqan3::sam_flag, (std::tuple_element_t<8, record_type>));
    EXPECT_SAME_TYPE(mate_t, (std::tuple_element_t<9, record_type>));
    EXPECT_SAME_TYPE(uint8_t, (std::tuple_element_t<10, record_type>));
    EXPECT_SAME_TYPE(std::vector<seqan3::cigar>, (std::tuple_element_t<11, record_type>));
    EXPECT_SAME_TYPE(seqan3::sam_tag_dictionary, (std::tuple_element_t<12, record_type>));

    EXPECT_EQ(std::tuple_size_v<record_type>, 13ul);

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
    EXPECT_DOUBLE_EQ(std::get<3>(r), 1);
    EXPECT_EQ(std::get<4>(r), alignment_t{});
    EXPECT_EQ(std::get<5>(r), "ref");
    EXPECT_EQ(std::get<6>(r), 0);
    EXPECT_EQ(std::get<7>(r), nullptr);
    EXPECT_EQ(std::get<8>(r), seqan3::sam_flag{41u});
    EXPECT_EQ(std::get<9>(r), (mate_t{0, 9, 300}));
    EXPECT_EQ(std::get<10>(r), 61u);
    EXPECT_RANGE_EQ(std::get<11>(r),
                    (std::vector<seqan3::cigar>{{1, 'S'_cigar_operation}, {1, 'M'_cigar_operation},
                                                {1, 'D'_cigar_operation}, {1, 'M'_cigar_operation},
                                                {1, 'I'_cigar_operation}}));
    EXPECT_EQ(std::get<12>(r), seqan3::sam_tag_dictionary{});
}

TEST_F(sam_record, get_by_type)
{
    record_type r{construct()};

    // EXPECT_EQ(std::get<std::string>(r), "MY ID"); // ambiguous
    EXPECT_RANGE_EQ(std::get<seqan3::dna5_vector>(r), "ACGT"_dna5);
    EXPECT_RANGE_EQ(std::get<std::vector<seqan3::phred42>>(r), "!##$"_phred42);
    EXPECT_DOUBLE_EQ(std::get<int32_t>(r), 1);
    EXPECT_EQ(std::get<alignment_t>(r), alignment_t{});
    // EXPECT_EQ(std::get<std::string>(r), "ref"); // ambiguous
    EXPECT_EQ(std::get<std::optional<int32_t>>(r), 0);
    EXPECT_EQ(std::get<seqan3::sam_file_header<std::vector<std::string>> *>(r), nullptr);
    EXPECT_EQ(std::get<seqan3::sam_flag>(r), seqan3::sam_flag{41u});
    EXPECT_EQ(std::get<mate_t>(r), (mate_t{0, 9, 300}));
    EXPECT_EQ(std::get<uint8_t>(r), 61u);
    EXPECT_RANGE_EQ(std::get<std::vector<seqan3::cigar>>(r),
                    (std::vector<seqan3::cigar>{{1, 'S'_cigar_operation}, {1, 'M'_cigar_operation},
                                                {1, 'D'_cigar_operation}, {1, 'M'_cigar_operation},
                                                {1, 'I'_cigar_operation}}));
    EXPECT_EQ(std::get<seqan3::sam_tag_dictionary>(r), seqan3::sam_tag_dictionary{});
}

TEST_F(sam_record, get_by_field)
{
    record_type r{construct()};

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
    EXPECT_EQ(seqan3::get<seqan3::field::id>(r), "MY ID");
    EXPECT_RANGE_EQ(seqan3::get<seqan3::field::seq>(r), "ACGT"_dna5);
    EXPECT_RANGE_EQ(seqan3::get<seqan3::field::qual>(r), "!##$"_phred42);
    EXPECT_DOUBLE_EQ(seqan3::get<seqan3::field::offset>(r), 1);
    EXPECT_EQ(seqan3::get<seqan3::field::alignment>(r), alignment_t{});
    EXPECT_EQ(seqan3::get<seqan3::field::ref_id>(r), "ref");
    EXPECT_EQ(seqan3::get<seqan3::field::ref_offset>(r), 0);
    EXPECT_EQ(seqan3::get<seqan3::field::header_ptr>(r), nullptr);
    EXPECT_EQ(seqan3::get<seqan3::field::flag>(r), seqan3::sam_flag{41u});
    EXPECT_EQ(seqan3::get<seqan3::field::mate>(r), (mate_t{0, 9, 300}));
    EXPECT_EQ(seqan3::get<seqan3::field::mapq>(r), 61u);
    EXPECT_RANGE_EQ(seqan3::get<seqan3::field::cigar>(r),
                    (std::vector<seqan3::cigar>{{1, 'S'_cigar_operation}, {1, 'M'_cigar_operation},
                                                {1, 'D'_cigar_operation}, {1, 'M'_cigar_operation},
                                                {1, 'I'_cigar_operation}}));
    EXPECT_EQ(seqan3::get<seqan3::field::tags>(r), seqan3::sam_tag_dictionary{});
#pragma GCC diagnostic pop
}

TEST_F(sam_record, get_by_member)
{
    record_type r{construct()};

    EXPECT_EQ(r.id(), "MY ID");
    EXPECT_RANGE_EQ(r.sequence(), "ACGT"_dna5);
    EXPECT_RANGE_EQ(r.base_qualities(), "!##$"_phred42);
    EXPECT_DOUBLE_EQ(r.sequence_position(), 1);
    EXPECT_EQ(r.alignment(), alignment_t{});
    EXPECT_EQ(r.reference_id(), "ref");
    EXPECT_EQ(r.reference_position(), 0);
    EXPECT_EQ(r.header_ptr(), nullptr);
    EXPECT_EQ(r.flag(), seqan3::sam_flag{41u});
    EXPECT_EQ(r.mate_reference_id(), 0);
    EXPECT_EQ(r.mate_position(), 9);
    EXPECT_EQ(r.template_length(), 300);
    EXPECT_EQ(r.mapping_quality(), 61u);
    EXPECT_RANGE_EQ(r.cigar_sequence(),
                    (std::vector<seqan3::cigar>{{1, 'S'_cigar_operation}, {1, 'M'_cigar_operation},
                                                {1, 'D'_cigar_operation}, {1, 'M'_cigar_operation},
                                                {1, 'I'_cigar_operation}}));
    EXPECT_EQ(r.tags(), seqan3::sam_tag_dictionary{});
}

TEST_F(sam_record, get_types)
{
    record_type r{construct()};

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
    EXPECT_SAME_TYPE(std::string &, decltype(seqan3::get<seqan3::field::id>(r)));
    EXPECT_SAME_TYPE(seqan3::dna5_vector &, decltype(seqan3::get<seqan3::field::seq>(r)));
    EXPECT_SAME_TYPE(std::vector<seqan3::phred42> &, decltype(seqan3::get<seqan3::field::qual>(r)));
    EXPECT_SAME_TYPE(int32_t &, decltype(seqan3::get<seqan3::field::offset>(r)));
    EXPECT_SAME_TYPE(alignment_t &, decltype(seqan3::get<seqan3::field::alignment>(r)));
    EXPECT_SAME_TYPE(std::string &, decltype(seqan3::get<seqan3::field::ref_id>(r)));
    EXPECT_SAME_TYPE(std::optional<int32_t> &, decltype(seqan3::get<seqan3::field::ref_offset>(r)));
    EXPECT_SAME_TYPE(seqan3::sam_file_header<std::vector<std::string>> * &,
                     decltype(seqan3::get<seqan3::field::header_ptr>(r)));
    EXPECT_SAME_TYPE(seqan3::sam_flag &, decltype(seqan3::get<seqan3::field::flag>(r)));
    EXPECT_SAME_TYPE(mate_t &, decltype(seqan3::get<seqan3::field::mate>(r)));
    EXPECT_SAME_TYPE(uint8_t &, decltype(seqan3::get<seqan3::field::mapq>(r)));
    EXPECT_SAME_TYPE(std::vector<seqan3::cigar> &, decltype(seqan3::get<seqan3::field::cigar>(r)));
    EXPECT_SAME_TYPE(seqan3::sam_tag_dictionary &, decltype(seqan3::get<seqan3::field::tags>(r)));

    EXPECT_SAME_TYPE(std::string const &, decltype(seqan3::get<seqan3::field::id>(std::as_const(r))));
    EXPECT_SAME_TYPE(seqan3::dna5_vector const &, decltype(seqan3::get<seqan3::field::seq>(std::as_const(r))));
    EXPECT_SAME_TYPE(std::vector<seqan3::phred42> const &,
                     decltype(seqan3::get<seqan3::field::qual>(std::as_const(r))));
    EXPECT_SAME_TYPE(int32_t const &, decltype(seqan3::get<seqan3::field::offset>(std::as_const(r))));
    EXPECT_SAME_TYPE(alignment_t const &, decltype(seqan3::get<seqan3::field::alignment>(std::as_const(r))));
    EXPECT_SAME_TYPE(std::string const &, decltype(seqan3::get<seqan3::field::ref_id>(std::as_const(r))));
    EXPECT_SAME_TYPE(std::optional<int32_t> const &,
                     decltype(seqan3::get<seqan3::field::ref_offset>(std::as_const(r))));
    EXPECT_SAME_TYPE(seqan3::sam_file_header<std::vector<std::string>> * const &,
                     decltype(seqan3::get<seqan3::field::header_ptr>(std::as_const(r))));
    EXPECT_SAME_TYPE(seqan3::sam_flag const &, decltype(seqan3::get<seqan3::field::flag>(std::as_const(r))));
    EXPECT_SAME_TYPE(mate_t const &, decltype(seqan3::get<seqan3::field::mate>(std::as_const(r))));
    EXPECT_SAME_TYPE(uint8_t const &, decltype(seqan3::get<seqan3::field::mapq>(std::as_const(r))));
    EXPECT_SAME_TYPE(std::vector<seqan3::cigar> const &, decltype(seqan3::get<seqan3::field::cigar>(std::as_const(r))));
    EXPECT_SAME_TYPE(seqan3::sam_tag_dictionary const &, decltype(seqan3::get<seqan3::field::tags>(std::as_const(r))));

    EXPECT_SAME_TYPE(std::string &&, decltype(seqan3::get<seqan3::field::id>(std::move(r))));
    EXPECT_SAME_TYPE(seqan3::dna5_vector &&, decltype(seqan3::get<seqan3::field::seq>(std::move(r))));
    EXPECT_SAME_TYPE(std::vector<seqan3::phred42> &&, decltype(seqan3::get<seqan3::field::qual>(std::move(r))));
    EXPECT_SAME_TYPE(int32_t &&, decltype(seqan3::get<seqan3::field::offset>(std::move(r))));
    EXPECT_SAME_TYPE(alignment_t &&, decltype(seqan3::get<seqan3::field::alignment>(std::move(r))));
    EXPECT_SAME_TYPE(std::string &&, decltype(seqan3::get<seqan3::field::ref_id>(std::move(r))));
    EXPECT_SAME_TYPE(std::optional<int32_t> &&, decltype(seqan3::get<seqan3::field::ref_offset>(std::move(r))));
    EXPECT_SAME_TYPE(seqan3::sam_file_header<std::vector<std::string>> * &&,
                     decltype(seqan3::get<seqan3::field::header_ptr>(std::move(r))));
    EXPECT_SAME_TYPE(seqan3::sam_flag &&, decltype(seqan3::get<seqan3::field::flag>(std::move(r))));
    EXPECT_SAME_TYPE(mate_t &&, decltype(seqan3::get<seqan3::field::mate>(std::move(r))));
    EXPECT_SAME_TYPE(uint8_t &&, decltype(seqan3::get<seqan3::field::mapq>(std::move(r))));
    EXPECT_SAME_TYPE(std::vector<seqan3::cigar> &&, decltype(seqan3::get<seqan3::field::cigar>(std::move(r))));
    EXPECT_SAME_TYPE(seqan3::sam_tag_dictionary &&, decltype(seqan3::get<seqan3::field::tags>(std::move(r))));

    EXPECT_SAME_TYPE(std::string const &&, decltype(seqan3::get<seqan3::field::id>(std::move(std::as_const(r)))));
    EXPECT_SAME_TYPE(seqan3::dna5_vector const &&,
                     decltype(seqan3::get<seqan3::field::seq>(std::move(std::as_const(r)))));
    EXPECT_SAME_TYPE(std::vector<seqan3::phred42> const &&,
                     decltype(seqan3::get<seqan3::field::qual>(std::move(std::as_const(r)))));
    EXPECT_SAME_TYPE(int32_t const &&, decltype(seqan3::get<seqan3::field::offset>(std::move(std::as_const(r)))));
    EXPECT_SAME_TYPE(alignment_t const &&,
                     decltype(seqan3::get<seqan3::field::alignment>(std::move(std::as_const(r)))));
    EXPECT_SAME_TYPE(std::string const &&, decltype(seqan3::get<seqan3::field::ref_id>(std::move(std::as_const(r)))));
    EXPECT_SAME_TYPE(std::optional<int32_t> const &&,
                     decltype(seqan3::get<seqan3::field::ref_offset>(std::move(std::as_const(r)))));
    EXPECT_SAME_TYPE(seqan3::sam_file_header<std::vector<std::string>> * const &&,
                     decltype(seqan3::get<seqan3::field::header_ptr>(std::move(std::as_const(r)))));
    EXPECT_SAME_TYPE(seqan3::sam_flag const &&,
                     decltype(seqan3::get<seqan3::field::flag>(std::move(std::as_const(r)))));
    EXPECT_SAME_TYPE(mate_t const &&, decltype(seqan3::get<seqan3::field::mate>(std::move(std::as_const(r)))));
    EXPECT_SAME_TYPE(uint8_t const &&, decltype(seqan3::get<seqan3::field::mapq>(std::move(std::as_const(r)))));
    EXPECT_SAME_TYPE(std::vector<seqan3::cigar> const &&,
                     decltype(seqan3::get<seqan3::field::cigar>(std::move(std::as_const(r)))));
    EXPECT_SAME_TYPE(seqan3::sam_tag_dictionary const &&,
                     decltype(seqan3::get<seqan3::field::tags>(std::move(std::as_const(r)))));
#pragma GCC diagnostic pop
}

TEST_F(sam_record, member_types)
{
    record_type r{construct()};

    EXPECT_SAME_TYPE(std::string &, decltype(r.id()));
    EXPECT_SAME_TYPE(seqan3::dna5_vector &, decltype(r.sequence()));
    EXPECT_SAME_TYPE(std::vector<seqan3::phred42> &, decltype(r.base_qualities()));
    EXPECT_SAME_TYPE(int32_t &, decltype(r.sequence_position()));
    EXPECT_SAME_TYPE(alignment_t &, decltype(r.alignment()));
    EXPECT_SAME_TYPE(std::string &, decltype(r.reference_id()));
    EXPECT_SAME_TYPE(std::optional<int32_t> &, decltype(r.reference_position()));
    EXPECT_SAME_TYPE(seqan3::sam_file_header<std::vector<std::string>> * &, decltype(r.header_ptr()));
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
    EXPECT_SAME_TYPE(int32_t const &, decltype(std::as_const(r.sequence_position())));
    EXPECT_SAME_TYPE(alignment_t const &, decltype(std::as_const(r.alignment())));
    EXPECT_SAME_TYPE(std::string const &, decltype(std::as_const(r.reference_id())));
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
    EXPECT_SAME_TYPE(int32_t &&, decltype(std::move(r.sequence_position())));
    EXPECT_SAME_TYPE(alignment_t &&, decltype(std::move(r.alignment())));
    EXPECT_SAME_TYPE(std::string &&, decltype(std::move(r.reference_id())));
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
    EXPECT_SAME_TYPE(int32_t const &&, decltype(std::move(std::as_const(r.sequence_position()))));
    EXPECT_SAME_TYPE(alignment_t const &&, decltype(std::move(std::as_const(r.alignment()))));
    EXPECT_SAME_TYPE(std::string const &&, decltype(std::move(std::as_const(r.reference_id()))));
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

// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <iostream>

#include <seqan3/alphabet/detail/debug_stream_alphabet.hpp>
#include <seqan3/core/debug_stream/byte.hpp>
#include <seqan3/core/debug_stream/optional.hpp>
#include <seqan3/core/debug_stream/tuple.hpp>
#include <seqan3/core/debug_stream/variant.hpp>
#include <seqan3/io/sam_file/input.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/test/fixture/io/sam_file/simple_three_verbose_reads_fixture.hpp>
#include <seqan3/test/pretty_printing.hpp>

using sam_file_seek_test_fixture = std::tuple<std::filesystem::path, std::vector<std::streampos>>;

struct sam_file_seek_test : public ::testing::TestWithParam<sam_file_seek_test_fixture>
{
    void SetUp() override
    {
        std::tie(sam_file_path, file_positions) = GetParam();

        sam_file_path = std::filesystem::path{CURRENT_SOURCE_DIR} / sam_file_path;
    }

    template <typename record_t, typename expected_record_t>
    static void expect_record_eq(record_t & record, expected_record_t & expected_record)
    {
        EXPECT_EQ(record.sequence(), expected_record.sequence());
        EXPECT_EQ(record.id(), expected_record.id());
        EXPECT_EQ(record.base_qualities(), expected_record.base_qualities());
        EXPECT_EQ(record.reference_id(), expected_record.reference_id());
        EXPECT_EQ(record.reference_position(), expected_record.reference_position());
        EXPECT_RANGE_EQ(record.cigar_sequence(), expected_record.cigar_sequence());
        EXPECT_EQ(record.flag(), expected_record.flag());
        EXPECT_EQ(record.mapping_quality(), expected_record.mapping_quality());
        EXPECT_EQ(record.mate_reference_id(), expected_record.mate_reference_id());
        EXPECT_EQ(record.mate_position(), expected_record.mate_position());
        EXPECT_EQ(record.template_length(), expected_record.template_length());
        EXPECT_EQ(record.tags(), expected_record.tags());
    }

    std::filesystem::path sam_file_path;
    std::vector<std::streampos> file_positions;
};

TEST_P(sam_file_seek_test, seek_to)
{
    seqan3::test::fixture::io::sam_file::simple_three_verbose_reads_fixture expected_file{};
    seqan3::sam_file_input fin{sam_file_path};

    ASSERT_GE(expected_file.records.size(), 3u);

    auto it = fin.begin();

    for (size_t i = 0u; i < expected_file.records.size(); ++it, ++i)
    {
        SCOPED_TRACE("sequential access");
        ASSERT_EQ(it.file_position(), file_positions[i]);
        expect_record_eq(*it, expected_file.records[i]);

        EXPECT_TRUE(it != fin.end());
    }
    EXPECT_TRUE(it == fin.end());

    for (size_t i : std::vector<size_t>{2u, 1u, 0u, 1u, 0u, 2u, 0u, 0u, 2u, 2u, 1u, 1u})
    {
        SCOPED_TRACE("random access");
        it.seek_to(file_positions[i]);
        expect_record_eq(*it, expected_file.records[i]);

        EXPECT_TRUE(it != fin.end());
    }
    EXPECT_TRUE(it != fin.end());

    for (size_t i = 1u; i < expected_file.records.size(); ++it, ++i)
    {
        SCOPED_TRACE("finish access sequentially");
        expect_record_eq(*it, expected_file.records[i]);

        EXPECT_TRUE(it != fin.end());
    }
    EXPECT_TRUE(it == fin.end());
}

INSTANTIATE_TEST_SUITE_P(bam_file,
                         sam_file_seek_test,
                         ::testing::Values(sam_file_seek_test_fixture{"simple_three_verbose_reads.bam",
                                                                      {4'915'200, 11'730'944, 23'134'208}}));

INSTANTIATE_TEST_SUITE_P(sam_file,
                         sam_file_seek_test,
                         ::testing::Values(sam_file_seek_test_fixture{"simple_three_verbose_reads.sam",
                                                                      {28, 135, 325}}));

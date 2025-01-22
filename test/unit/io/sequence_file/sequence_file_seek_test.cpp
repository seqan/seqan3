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
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/test/fixture/io/sequence_file/standard_fixture.hpp>
#include <seqan3/test/pretty_printing.hpp>

using sequence_file_seek_test_fixture = std::tuple<std::filesystem::path /*sequence_file_path*/,
                                                   bool /*has_base_qualities*/,
                                                   std::vector<std::streampos> /*file_positions*/>;

struct sequence_file_seek_test : public ::testing::TestWithParam<sequence_file_seek_test_fixture>
{
    void SetUp() override
    {
        std::tie(sequence_file_path, has_base_qualities, file_positions) = GetParam();

        sequence_file_path = std::filesystem::path{CURRENT_SOURCE_DIR} / sequence_file_path;
    }

    template <typename record_t, typename expected_record_t>
    void expect_record_eq(record_t & record, expected_record_t & expected_record)
    {
        EXPECT_EQ(record.sequence(), expected_record.sequence());
        EXPECT_EQ(record.id(), expected_record.id());
        if (has_base_qualities)
        {
            EXPECT_EQ(record.base_qualities(), expected_record.base_qualities());
        }
    }

    std::filesystem::path sequence_file_path;
    bool has_base_qualities;
    std::vector<std::streampos> file_positions;
};

TEST_P(sequence_file_seek_test, seek_to)
{
    seqan3::test::fixture::io::sequence_file::standard_fixture expected_file{};
    seqan3::sequence_file_input fin{sequence_file_path};

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

sequence_file_seek_test_fixture fasta_file_fixture{"standard.fasta", false, {0, 25, 114}};
INSTANTIATE_TEST_SUITE_P(fasta_file, sequence_file_seek_test, ::testing::Values(fasta_file_fixture));

sequence_file_seek_test_fixture fastq_file_fixture{"standard.fastq", true, {0, 45, 218}};
INSTANTIATE_TEST_SUITE_P(fastq_file, sequence_file_seek_test, ::testing::Values(fastq_file_fixture));

sequence_file_seek_test_fixture sam_file_fixture{"standard.sam", true, {49, 107, 293}};
INSTANTIATE_TEST_SUITE_P(sam_file, sequence_file_seek_test, ::testing::Values(sam_file_fixture));

sequence_file_seek_test_fixture embl_file_fixture{"standard.embl", false, {0, 108, 283}};
INSTANTIATE_TEST_SUITE_P(embl_file, sequence_file_seek_test, ::testing::Values(embl_file_fixture));

sequence_file_seek_test_fixture genbank_file_fixture{"standard.genbank", false, {0, 561, 802}};
INSTANTIATE_TEST_SUITE_P(genbank_file, sequence_file_seek_test, ::testing::Values(genbank_file_fixture));

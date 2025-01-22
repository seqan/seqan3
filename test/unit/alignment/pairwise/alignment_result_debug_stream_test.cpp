// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <utility>
#include <vector>

#include <seqan3/alignment/aligned_sequence/debug_stream_alignment.hpp>
#include <seqan3/alignment/pairwise/alignment_result.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream/tuple.hpp>

TEST(alignment_result_test, debug_streamable)
{
    using coordinate_t = std::pair<int, int>;
    using aligned_seq_type = std::vector<seqan3::gapped<seqan3::dna4>>;

    using seqan3::operator""_dna4;

    size_t id = 3;
    int score = -15;
    coordinate_t begin_coordinate{4, 6};
    coordinate_t end_coordinate{23, 35};

    aligned_seq_type gapped_seq1{'A'_dna4, 'T'_dna4, seqan3::gap{}, 'A'_dna4};
    aligned_seq_type gapped_seq2{'A'_dna4, 'T'_dna4, 'C'_dna4, seqan3::gap{}};
    std::pair<aligned_seq_type, aligned_seq_type> alignment{gapped_seq1, gapped_seq2};

    std::ostringstream ostream{};
    seqan3::debug_stream_type debug_stream{ostream};

    { // Print id and score
        seqan3::detail::alignment_result_value_type result_value{id, id, score};
        seqan3::alignment_result result{result_value};
        debug_stream << result;

        EXPECT_EQ(ostream.str(), "{sequence1 id: 3, sequence2 id: 3, score: -15}");
    }

    { // Print id and score and back coordinate
        ostream.str("");
        seqan3::detail::alignment_result_value_type result_value{id, id, score, end_coordinate};
        seqan3::alignment_result result{result_value};

        EXPECT_TRUE((seqan3::detail::is_type_specialisation_of_v<decltype(result), seqan3::alignment_result>));
        debug_stream << result;

        EXPECT_EQ(ostream.str(), "{sequence1 id: 3, sequence2 id: 3, score: -15, end: (23,35)}");
    }

    { // Print id and score and back and front coordinate
        ostream.str("");
        seqan3::detail::alignment_result_value_type result_value{id, id, score, end_coordinate, begin_coordinate};
        seqan3::alignment_result result{result_value};
        debug_stream << result;

        EXPECT_EQ(ostream.str(), "{sequence1 id: 3, sequence2 id: 3, score: -15, begin: (4,6), end: (23,35)}");
    }

    { // Print id and score and back and front coordinate and alignment
        ostream.str("");
        seqan3::detail::alignment_result_value_type result_value{id,
                                                                 id,
                                                                 score,
                                                                 end_coordinate,
                                                                 begin_coordinate,
                                                                 alignment};
        seqan3::alignment_result result{result_value};
        debug_stream << result;

        EXPECT_EQ(ostream.str(),
                  "{sequence1 id: 3, sequence2 id: 3, score: -15, begin: (4,6), end: (23,35), \n"
                  "alignment:\n"
                  "      0     \n"
                  "        AT-A\n"
                  "        ||  \n"
                  "        ATC-\n"
                  "}");
    }
}

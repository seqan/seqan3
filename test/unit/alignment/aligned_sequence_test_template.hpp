// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <iterator>
#include <string>

#include <seqan3/alignment/aligned_sequence/aligned_sequence_concept.hpp>
#include <seqan3/alignment/cigar_conversion/cigar_from_alignment.hpp>
#include <seqan3/alphabet/detail/debug_stream_alphabet.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream/detail/to_string.hpp>
#include <seqan3/core/debug_stream/range.hpp>
#include <seqan3/io/sam_file/detail/cigar.hpp>
#include <seqan3/test/expect_range_eq.hpp>

using seqan3::operator""_dna4;

template <typename t>
class aligned_sequence : public ::testing::Test
{};

seqan3::dna4_vector const seq = "ACTA"_dna4;

TYPED_TEST_SUITE_P(aligned_sequence);

TYPED_TEST_P(aligned_sequence, fulfills_concept)
{
    EXPECT_TRUE((seqan3::aligned_sequence<TypeParam>));
    EXPECT_FALSE((seqan3::aligned_sequence<std::vector<seqan3::dna4>>));
}

TYPED_TEST_P(aligned_sequence, assign_unaligned_sequence)
{
    using unaligned_seq_type = std::remove_cvref_t<seqan3::detail::unaligned_seq_t<TypeParam>>;
    unaligned_seq_type unaligned_seq{};

    if constexpr (seqan3::sequence_container<unaligned_seq_type>)
    {
        unaligned_seq.resize(seq.size());
        std::copy(seq.begin(), seq.end(), std::ranges::begin(unaligned_seq));
    }
    else // type is view, happens for gap_decorator tests
        unaligned_seq = unaligned_seq_type{seq};

    TypeParam aligned_seq{};

    assign_unaligned(aligned_seq, unaligned_seq);

    EXPECT_EQ(aligned_seq.size(), unaligned_seq.size());
    EXPECT_RANGE_EQ(aligned_seq, unaligned_seq);
}

TYPED_TEST_P(aligned_sequence, assign_empty_unaligned_sequence)
{
    using unaligned_seq_type = std::remove_cvref_t<seqan3::detail::unaligned_seq_t<TypeParam>>;
    unaligned_seq_type unaligned_seq{};
    TypeParam aligned_seq{};

    assign_unaligned(aligned_seq, unaligned_seq);

    EXPECT_EQ(aligned_seq.size(), unaligned_seq.size());
    EXPECT_EQ(aligned_seq.size(), 0u);
}

TYPED_TEST_P(aligned_sequence, insert_one_gap)
{
    TypeParam aligned_seq;
    TestFixture::initialise_typed_test_container(aligned_seq, seq);

    EXPECT_EQ(aligned_seq.size(), 4u);

    auto it = seqan3::insert_gap(aligned_seq, std::ranges::next(std::ranges::begin(aligned_seq), 1));
    EXPECT_EQ(*it, seqan3::gap{});
    EXPECT_EQ(aligned_seq[1], seqan3::gap{});
    EXPECT_EQ(aligned_seq.size(), 5u);
    EXPECT_EQ(seqan3::detail::to_string(aligned_seq), "A-CTA");

    it = seqan3::insert_gap(aligned_seq, std::ranges::next(std::ranges::begin(aligned_seq), 1));
    EXPECT_EQ(*it, seqan3::gap{});
    EXPECT_EQ(aligned_seq[1], seqan3::gap{});
    EXPECT_EQ(aligned_seq[2], seqan3::gap{});
    EXPECT_EQ(aligned_seq.size(), 6u);
    EXPECT_EQ(seqan3::detail::to_string(aligned_seq), "A--CTA");
}

TYPED_TEST_P(aligned_sequence, insert_multiple_gaps)
{
    TypeParam aligned_seq;
    TestFixture::initialise_typed_test_container(aligned_seq, seq);
    EXPECT_EQ(aligned_seq.size(), 4u);

    auto it = seqan3::insert_gap(aligned_seq, std::ranges::next(begin(aligned_seq), 1), 2);
    EXPECT_EQ(*it, seqan3::gap{});
    EXPECT_EQ(*++it, seqan3::gap{});
    EXPECT_EQ(aligned_seq[1], seqan3::gap{});
    EXPECT_EQ(aligned_seq[2], seqan3::gap{});
    EXPECT_EQ(aligned_seq.size(), 6u);

    // insert a gap within another gap
    seqan3::insert_gap(aligned_seq, std::ranges::next(begin(aligned_seq), 2), 4);
    EXPECT_EQ(seqan3::detail::to_string(aligned_seq), "A------CTA");

    // insert at begin
    seqan3::insert_gap(aligned_seq, std::ranges::begin(aligned_seq), 2);
    EXPECT_EQ(seqan3::detail::to_string(aligned_seq), "--A------CTA");

    // insert at end
    seqan3::insert_gap(aligned_seq, std::ranges::end(aligned_seq), 2);
    EXPECT_EQ(seqan3::detail::to_string(aligned_seq), "--A------CTA--");
}

TYPED_TEST_P(aligned_sequence, insert_zero_gaps)
{
    TypeParam aligned_seq;
    TestFixture::initialise_typed_test_container(aligned_seq, seq);
    EXPECT_EQ(aligned_seq.size(), 4u);

    auto it = seqan3::insert_gap(aligned_seq, std::ranges::next(begin(aligned_seq), 1), 0);
    typename TypeParam::value_type val{'C'_dna4};
    EXPECT_EQ(*it, val);
    EXPECT_EQ(seqan3::detail::to_string(aligned_seq), "ACTA");
}

TYPED_TEST_P(aligned_sequence, erase_one_gap)
{
    // 1) Removing an actual gap
    TypeParam aligned_seq;
    TestFixture::initialise_typed_test_container(aligned_seq, seq);
    EXPECT_EQ(aligned_seq.size(), 4u);
    EXPECT_EQ(seqan3::detail::to_string(aligned_seq), "ACTA");

    insert_gap(aligned_seq, std::ranges::next(begin(aligned_seq), 1));
    EXPECT_EQ(aligned_seq.size(), 5u);
    EXPECT_EQ(seqan3::detail::to_string(aligned_seq), "A-CTA");

    typename TypeParam::value_type val{'C'_dna4};
    auto it = seqan3::erase_gap(aligned_seq, std::ranges::next(begin(aligned_seq), 1));
    EXPECT_EQ(*it, val);
    EXPECT_EQ(aligned_seq.size(), 4u);
    EXPECT_EQ(seqan3::detail::to_string(aligned_seq), "ACTA");

    // 2) Removing a non-gap
    EXPECT_THROW(seqan3::erase_gap(aligned_seq, std::ranges::next(std::ranges::begin(aligned_seq), 2)),
                 seqan3::gap_erase_failure);
    EXPECT_EQ(aligned_seq.size(), 4u); // no change
    EXPECT_EQ(seqan3::detail::to_string(aligned_seq), "ACTA");
}

TYPED_TEST_P(aligned_sequence, erase_multiple_gaps)
{
    TypeParam aligned_seq;

    // 1) Removing a gap of length > 1
    TestFixture::initialise_typed_test_container(aligned_seq, seq);
    EXPECT_EQ(aligned_seq.size(), 4u);
    EXPECT_EQ(seqan3::detail::to_string(aligned_seq), "ACTA");

    seqan3::insert_gap(aligned_seq, std::ranges::next(std::ranges::begin(aligned_seq), 1), 2);
    EXPECT_EQ(aligned_seq.size(), 6u);
    EXPECT_EQ(seqan3::detail::to_string(aligned_seq), "A--CTA");

    typename TypeParam::value_type val{'C'_dna4};
    auto it = seqan3::erase_gap(aligned_seq,
                                std::ranges::next(std::ranges::begin(aligned_seq), 1),
                                std::ranges::next(std::ranges::begin(aligned_seq), 3));
    EXPECT_EQ(aligned_seq.size(), 4u);
    EXPECT_EQ(*it, val);
    EXPECT_EQ(seqan3::detail::to_string(aligned_seq), "ACTA");

    // 2) Removing a non-gap
    TestFixture::initialise_typed_test_container(aligned_seq, seq); // reset

    EXPECT_THROW(seqan3::erase_gap(aligned_seq,
                                   std::ranges::next(std::ranges::begin(aligned_seq), 1),
                                   std::ranges::next(std::ranges::begin(aligned_seq), 3)),
                 seqan3::gap_erase_failure);
    EXPECT_EQ(aligned_seq.size(), seq.size()); // no change
    EXPECT_EQ(seqan3::detail::to_string(aligned_seq), "ACTA");

    // 3) Remove a gap of length 1 within a gap of length 5
    TestFixture::initialise_typed_test_container(aligned_seq, seq); // reset
    insert_gap(aligned_seq, std::ranges::next(std::ranges::begin(aligned_seq), 1), 5);
    EXPECT_EQ(aligned_seq.size(), seq.size() + 5);
    EXPECT_EQ(seqan3::detail::to_string(aligned_seq), "A-----CTA");

    it = erase_gap(aligned_seq, std::ranges::next(std::ranges::begin(aligned_seq), 3)); // erase one in the middle of 5
    EXPECT_EQ(aligned_seq.size(), seq.size() + 4);
    EXPECT_EQ(aligned_seq[5], val);
    EXPECT_EQ(*it, seqan3::gap{});
    EXPECT_EQ(seqan3::detail::to_string(aligned_seq), "A----CTA");

    // 4) Remove gaps two times
    TestFixture::initialise_typed_test_container(aligned_seq, seq); // reset
    insert_gap(aligned_seq, std::ranges::next(std::ranges::begin(aligned_seq), 3), 4);
    insert_gap(aligned_seq, std::ranges::next(std::ranges::begin(aligned_seq), 1), 5);
    EXPECT_EQ(aligned_seq.size(), seq.size() + 9);
    EXPECT_EQ(seqan3::detail::to_string(aligned_seq), "A-----CT----A");

    seqan3::erase_gap(aligned_seq,
                      std::ranges::next(std::ranges::begin(aligned_seq), 2),
                      std::ranges::next(std::ranges::begin(aligned_seq), 4));
    seqan3::erase_gap(aligned_seq,
                      std::ranges::next(std::ranges::begin(aligned_seq), 6),
                      std::ranges::next(std::ranges::begin(aligned_seq), 10));
    EXPECT_EQ(seqan3::detail::to_string(aligned_seq), "A---CTA");

    // 5) Removing too much from a gap
    EXPECT_THROW(seqan3::erase_gap(aligned_seq,
                                   std::ranges::next(std::ranges::begin(aligned_seq), 2),
                                   std::ranges::next(std::ranges::begin(aligned_seq), 5)),
                 seqan3::gap_erase_failure);
    EXPECT_EQ(aligned_seq.size(), 7u); // no change
    EXPECT_EQ(seqan3::detail::to_string(aligned_seq), "A---CTA");
}

TYPED_TEST_P(aligned_sequence, insert_erase_on_empty_sequence)
{
    using unaligned_seq_type = std::remove_cvref_t<seqan3::detail::unaligned_seq_t<TypeParam>>;
    unaligned_seq_type unaligned{};
    TypeParam aligned_seq{};

    assign_unaligned(aligned_seq, unaligned);

    auto it = insert_gap(aligned_seq, std::ranges::begin(aligned_seq));
    EXPECT_EQ(*it, seqan3::gap{});
    EXPECT_EQ(aligned_seq.size(), 1u);
    EXPECT_EQ(seqan3::detail::to_string(aligned_seq), "-");

    it = seqan3::insert_gap(aligned_seq, std::ranges::end(aligned_seq), 3);
    EXPECT_EQ(*it, seqan3::gap{});
    EXPECT_EQ(aligned_seq.size(), 4u);
    EXPECT_EQ(seqan3::detail::to_string(aligned_seq), "----");

    it = seqan3::insert_gap(aligned_seq, std::ranges::end(aligned_seq), 0);
    EXPECT_EQ(aligned_seq.size(), 4u);
    EXPECT_EQ(seqan3::detail::to_string(aligned_seq), "----");

    it = seqan3::erase_gap(aligned_seq, std::ranges::begin(aligned_seq));
    EXPECT_EQ(aligned_seq.size(), 3u);
    EXPECT_EQ(seqan3::detail::to_string(aligned_seq), "---");

    it = seqan3::erase_gap(aligned_seq, std::ranges::begin(aligned_seq), std::ranges::end(aligned_seq));
    EXPECT_EQ(aligned_seq.size(), 0u);
    EXPECT_EQ(seqan3::detail::to_string(aligned_seq), "");
}

TYPED_TEST_P(aligned_sequence, cigar_string)
{
    { // default_parameters
        auto seq_ref = "ACGTGATCTG"_dna4;
        auto seq_read = "ACGTCGTAGTG"_dna4;
        TypeParam ref;
        TypeParam read;
        TestFixture::initialise_typed_test_container(ref, seq_ref);
        TestFixture::initialise_typed_test_container(read, seq_read);
        seqan3::insert_gap(ref, std::ranges::next(begin(ref), 7), 2);
        seqan3::insert_gap(read, std::ranges::next(begin(read), 4), 1);

        auto cigar = seqan3::cigar_from_alignment(std::make_pair(ref, read));
        EXPECT_EQ(seqan3::detail::get_cigar_string(cigar), "4M1D2M2I3M");

        TypeParam ref2;
        TypeParam read2;
        TestFixture::initialise_typed_test_container(ref2, seq_ref);
        TestFixture::initialise_typed_test_container(read2, seq_read);
        insert_gap(ref2, std::ranges::next(begin(ref2), 10), 2);
        insert_gap(ref2, std::ranges::next(begin(ref2), 7), 2);
        insert_gap(ref2, std::ranges::begin(ref2), 3);
        insert_gap(read2, std::ranges::next(begin(read2), 11), 4);
        insert_gap(read2, std::ranges::next(begin(read2), 4), 1);
        insert_gap(read2, std::ranges::begin(read2), 1);

        auto cigar2 = seqan3::cigar_from_alignment(std::make_pair(ref2, read2));
        EXPECT_EQ(seqan3::detail::get_cigar_string(cigar2), "1P2I2M1D4M2I1M2D2P");
    }
    {
        // with soft clipping
        auto seq_ref = "ACGTGATCTG"_dna4;
        auto seq_read = "ACGTCGTAGTG"_dna4;
        TypeParam ref;
        TypeParam read;
        TestFixture::initialise_typed_test_container(ref, seq_ref);
        TestFixture::initialise_typed_test_container(read, seq_read);
        insert_gap(ref, std::ranges::next(std::ranges::begin(ref), 7), 2);
        insert_gap(read, std::ranges::next(std::ranges::begin(read), 4), 1);

        auto cigar = seqan3::cigar_from_alignment(std::make_pair(ref, read), {.soft_front = 5, .soft_back = 60});
        EXPECT_EQ(seqan3::detail::get_cigar_string(cigar), "5S4M1D2M2I3M60S");

        // gaps at the end
        TypeParam ref2;
        TypeParam read2;
        TestFixture::initialise_typed_test_container(ref2, seq_ref);
        TestFixture::initialise_typed_test_container(read2, seq_read);
        seqan3::insert_gap(ref2, std::ranges::next(begin(ref2), 10), 2);
        seqan3::insert_gap(ref2, std::ranges::next(begin(ref2), 7), 2);
        seqan3::insert_gap(ref2, std::ranges::begin(ref2), 3);
        seqan3::insert_gap(read2, std::ranges::next(begin(read2), 11), 4);
        seqan3::insert_gap(read2, std::ranges::next(begin(read2), 4), 1);
        seqan3::insert_gap(read2, std::ranges::begin(read2), 1);

        auto cigar2 = seqan3::cigar_from_alignment(std::make_pair(ref2, read2), {.soft_front = 3, .soft_back = 5});
        EXPECT_EQ(seqan3::detail::get_cigar_string(cigar2), "3S1P2I2M1D4M2I1M2D2P5S");
    }
    {
        // no gaps at the end
        auto seq_ref = "ACGTGATCAG"_dna4;
        auto seq_read = "ACGTCGTACTG"_dna4;
        TypeParam ref;
        TypeParam read;
        TestFixture::initialise_typed_test_container(ref, seq_ref);
        TestFixture::initialise_typed_test_container(read, seq_read);
        seqan3::insert_gap(ref, std::ranges::next(std::ranges::begin(ref), 7), 2);
        seqan3::insert_gap(read, std::ranges::next(std::ranges::begin(read), 4), 1);

        auto cigar = seqan3::cigar_from_alignment(std::make_pair(ref, read), {}, true);
        EXPECT_EQ(seqan3::detail::get_cigar_string(cigar), "4=1D2X2I1=1X1=");
        auto cigar2 = seqan3::cigar_from_alignment(std::make_pair(ref, read), {.soft_front = 5, .soft_back = 60}, true);
        EXPECT_EQ(seqan3::detail::get_cigar_string(cigar2), "5S4=1D2X2I1=1X1=60S");
    }
}

REGISTER_TYPED_TEST_SUITE_P(aligned_sequence,
                            fulfills_concept,
                            assign_unaligned_sequence,
                            assign_empty_unaligned_sequence,
                            insert_erase_on_empty_sequence,
                            insert_one_gap,
                            insert_multiple_gaps,
                            insert_zero_gaps,
                            erase_one_gap,
                            erase_multiple_gaps,
                            cigar_string);

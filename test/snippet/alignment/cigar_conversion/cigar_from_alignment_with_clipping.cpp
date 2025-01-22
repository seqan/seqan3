// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alignment/aligned_sequence/aligned_sequence_concept.hpp>
#include <seqan3/alignment/cigar_conversion/cigar_from_alignment.hpp>
#include <seqan3/alignment/decorator/gap_decorator.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/views/slice.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::dna5_vector reference = "ATGGCGTAGAGCTTCCCCCCCCCCCCCCCCC"_dna5;
    seqan3::dna5_vector read = "ATGCCCCGTTGCTT"_dna5; // length 14

    // Let's say, we want to ignore the last 2 bases of the query because the quality is low.
    // We thus only align the first 12 bases, the last two will be soft-clipped bases in the CIGAR string.
    seqan3::gap_decorator aligned_reference{reference | seqan3::views::slice(0, 12)};
    seqan3::gap_decorator aligned_query{read | seqan3::views::slice(0, 12)};
    // insert gaps
    seqan3::insert_gap(aligned_reference, aligned_reference.begin() + 4, 2);
    seqan3::insert_gap(aligned_query, aligned_query.begin() + 11, 2);

    auto cigar_sequence =
        seqan3::cigar_from_alignment(std::tie(aligned_reference, aligned_query),
                                     {.hard_front = 1, .hard_back = 0, .soft_front = 0, .soft_back = 2});

    seqan3::debug_stream << cigar_sequence << '\n'; // prints [1H,4M,2I,5M,2D,1M,2S]
}

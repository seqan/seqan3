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

    // Align the full query against the first 14 bases of the reference.
    seqan3::gap_decorator aligned_reference{reference | seqan3::views::slice(0, 14)};
    seqan3::gap_decorator aligned_read{read};
    // Insert gaps to represent the alignment:
    seqan3::insert_gap(aligned_read, aligned_read.begin() + 11, 2);
    seqan3::insert_gap(aligned_reference, aligned_reference.begin() + 4, 2);

    seqan3::debug_stream << aligned_reference << '\n' << aligned_read << '\n';
    // prints:
    // ATGG--CGTAGAGCTT
    // ATGCCCCGTTG--CTT

    auto cigar_sequence = seqan3::cigar_from_alignment(std::tie(aligned_reference, aligned_read));

    seqan3::debug_stream << cigar_sequence << '\n'; // prints [4M,2I,5M,2D,3M]
}

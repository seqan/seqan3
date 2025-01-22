// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alignment/cigar_conversion/alignment_from_cigar.hpp>
#include <seqan3/alphabet/cigar/cigar.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/core/debug_stream.hpp>

using namespace seqan3::literals;

int main()
{
    // CIGAR string = 2M1D2M
    std::vector<seqan3::cigar> cigar_vector{{2, 'M'_cigar_operation},
                                            {1, 'D'_cigar_operation},
                                            {2, 'M'_cigar_operation}};

    uint32_t reference_start_position{0}; // The read is aligned at the start of the reference.
    seqan3::dna5_vector reference = "ACTGATCGAGAGGATCTAGAGGAGATCGTAGGAC"_dna5;
    seqan3::dna5_vector query = "ACGA"_dna5;

    auto alignment = alignment_from_cigar(cigar_vector, reference, reference_start_position, query);

    seqan3::debug_stream << alignment << '\n'; // prints (ACTGA,AC-GA)
}

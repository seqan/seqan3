// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alignment/cigar_conversion/alignment_from_cigar.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sam_file/all.hpp>

using namespace seqan3::literals;

auto sam_file_raw = R"(@HD	VN:1.6
@SQ	SN:ref	LN:34
read1	41	ref	1	61	1S1M1D1M1I	ref	10	300	ACGT	!##$	AS:i:2	NM:i:7
read2	42	ref	2	62	1H7M1D1M1S2H	ref	10	300	AGGCTGNAG	!##$&'()*	xy:B:S,3,4,5
read3	43	ref	3	63	1S1M1P1M1I1M1I1D1M1S	ref	10	300	GGAGTATA	!!*+,-./
)";

int main()
{
    // The reference sequence might be read from a different file.
    seqan3::dna5_vector reference = "ACTGATCGAGAGGATCTAGAGGAGATCGTAGGAC"_dna5;

    seqan3::sam_file_input fin{std::istringstream{sam_file_raw}, seqan3::format_sam{}};
    // You will probably read it from a file, e.g., like this:
    // seqan3::sam_file_input fin{"test.sam"};

    for (auto && rec : fin)
    {
        auto alignment =
            alignment_from_cigar(rec.cigar_sequence(), reference, rec.reference_position().value(), rec.sequence());
        seqan3::debug_stream << alignment << '\n';
    }

    // prints:
    // (ACT-,C-GT)
    // (CTGATCGAG,AGGCTGN-A)
    // (T-G-A-TC,G-AGTA-T)
}

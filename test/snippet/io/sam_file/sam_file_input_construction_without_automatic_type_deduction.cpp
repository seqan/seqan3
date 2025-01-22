// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <sstream>

#include <seqan3/io/sam_file/input.hpp>
#include <seqan3/utility/type_list/type_list.hpp>

auto input = R"(@HD	VN:1.6	SO:coordinate
r001	99	ref	7	30	8M2I4M1D3M	=	37	39	TTAGATAAAGGATACTG	*)";

int main()
{
    // The default types; you can adjust this list if you don't want to read all this data.
    using default_fields = seqan3::fields<seqan3::field::seq,
                                          seqan3::field::id,
                                          seqan3::field::ref_id,
                                          seqan3::field::ref_offset,
                                          seqan3::field::cigar,
                                          seqan3::field::mapq,
                                          seqan3::field::qual,
                                          seqan3::field::flag,
                                          seqan3::field::mate,
                                          seqan3::field::tags,
                                          seqan3::field::header_ptr>;

    // The expected format:
    using sam_file_input_t = seqan3::sam_file_input<seqan3::sam_file_input_default_traits<>,
                                                    default_fields,
                                                    // Which formats are allowed:
                                                    seqan3::type_list<seqan3::format_sam>>;

    sam_file_input_t fin{std::istringstream{input}, seqan3::format_sam{}};
}

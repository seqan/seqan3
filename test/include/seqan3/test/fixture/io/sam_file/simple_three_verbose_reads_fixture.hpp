// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <optional>
#include <tuple>

#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/io/sam_file/header.hpp>
#include <seqan3/io/sam_file/record.hpp>
#include <seqan3/io/sam_file/sam_flag.hpp>
#include <seqan3/io/sam_file/sam_tag_dictionary.hpp>
#include <seqan3/test/literal/cigar_literal.hpp>
#include <seqan3/test/literal/gapped_dna5_literal.hpp>

namespace seqan3::test::fixture::io::sam_file
{

// "parsed" data from test/unit/io/sam_file/simple_three_verbose_reads.sam / bam
struct simple_three_verbose_reads_fixture
{
    using alignment_t = std::pair<std::vector<seqan3::gapped<seqan3::dna5>>, std::vector<seqan3::gapped<seqan3::dna5>>>;

    using mate_t = std::tuple<std::optional<int32_t>, std::optional<int32_t>, int32_t>;

    using header_type = seqan3::sam_file_header<std::vector<std::string>>;

    using types = seqan3::type_list<std::string,                  // seqan3::field::id,
                                    seqan3::dna5_vector,          // seqan3::field::seq,
                                    std::vector<seqan3::phred42>, // seqan3::field::qual,
                                    alignment_t,                  // seqan3::field::alignment,
                                    std::optional<int32_t>,       // seqan3::field::ref_id,
                                    std::optional<int32_t>,       // seqan3::field::ref_offset,
                                    header_type *,                // seqan3::field::header_ptr,
                                    seqan3::sam_flag,             // seqan3::field::flag,
                                    mate_t,                       // seqan3::field::mate,
                                    uint8_t,                      // seqan3::field::mapq,
                                    std::vector<seqan3::cigar>,   // seqan3::field::cigar,
                                    seqan3::sam_tag_dictionary>;  // seqan3::field::tags

    using types_as_ids = seqan3::fields<seqan3::field::id,
                                        seqan3::field::seq,
                                        seqan3::field::qual,
                                        seqan3::field::alignment,
                                        seqan3::field::ref_id,
                                        seqan3::field::ref_offset,
                                        seqan3::field::header_ptr,
                                        seqan3::field::flag,
                                        seqan3::field::mate,
                                        seqan3::field::mapq,
                                        seqan3::field::cigar,
                                        seqan3::field::tags>;

    using record_type = seqan3::sam_record<types, types_as_ids>;

    std::vector<std::string> reference_ids{"ref"};
    std::vector<seqan3::dna5_vector> reference_sequences{"ACTGATCGAGAGGATCTAGAGGAGATCGTAGGAC"_dna5};

    header_type header{[](auto && reference_ids) -> header_type
                       {
                           header_type header{reference_ids};
                           header.ref_id_info.emplace_back(reference_ids.size(), "");
                           header.ref_dict[header.ref_ids()[0]] =
                               0; // set up header which is otherwise done on file level
                           return header;
                       }(reference_ids)};

    record_type record1{/*.id =*/"read1",
                        /*.sequence =*/"ACGT"_dna5,
                        /*.base_qualities =*/"!##$"_phred42,
                        /*.alignment =*/alignment_t{"ACT-"_gapped_dna5, "C-GT"_gapped_dna5},
                        /*.reference_id =*/0, // "ref"
                        /*.reference_position =*/0,
                        /*.header_ptr =*/&header,
                        /*.flag =*/seqan3::sam_flag{41u},
                        {/*.mate_reference_id =*/0, /*.mate_position =*/9, /*.template_length =*/300},
                        /*.mapping_quality =*/61u,
                        /*.cigar_sequence =*/"1S1M1D1M1I"_cigar,
                        /*.tags =*/
                        seqan3::sam_tag_dictionary{[]()
                                                   {
                                                       seqan3::sam_tag_dictionary tags{};
                                                       tags["NM"_tag] = -7;
                                                       tags["AS"_tag] = 2;
                                                       tags["CC"_tag] = 300;
                                                       tags["cc"_tag] = -300;
                                                       tags["aa"_tag] = 'c';
                                                       tags["ff"_tag] = 3.1f;
                                                       tags["zz"_tag] = "str";
                                                       return tags;
                                                   }()}};

    record_type record2{/*.id =*/"read2",
                        /*.sequence =*/"AGGCTGNAG"_dna5,
                        /*.base_qualities =*/"!##$&'()*"_phred42,
                        /*.alignment =*/alignment_t{"CTGATCGAG"_gapped_dna5, "AGGCTGN-A"_gapped_dna5},
                        /*.reference_id =*/0, // "ref"
                        /*.reference_position =*/1,
                        /*.header_ptr =*/&header,
                        /*.flag =*/seqan3::sam_flag{42u},
                        {/*.mate_reference_id =*/0, /*.mate_position =*/9, /*.template_length =*/300},
                        /*.mapping_quality =*/62u,
                        /*.cigar_sequence =*/"1H7M1D1M1S"_cigar,
                        /*.tags =*/
                        seqan3::sam_tag_dictionary{[]()
                                                   {
                                                       seqan3::sam_tag_dictionary tags{};
                                                       tags["bc"_tag] = std::vector<int8_t>{-3};
                                                       tags["bC"_tag] = std::vector<uint8_t>{3u, 200u};
                                                       tags["bs"_tag] = std::vector<int16_t>{-3, 200, -300};
                                                       tags["bS"_tag] = std::vector<uint16_t>{300u, 40u, 500u};
                                                       tags["bi"_tag] = std::vector<int32_t>{-3, 200, -66000};
                                                       tags["bI"_tag] = std::vector<uint32_t>{294'967'296u};
                                                       tags["bf"_tag] = std::vector<float>{3.5f, 0.1f, 43.8f};
                                                       tags["bH"_tag] = std::vector<std::byte>{std::byte{0x1A},
                                                                                               std::byte{0xE3},
                                                                                               std::byte{0x01}};
                                                       return tags;
                                                   }()}};

    record_type record3{/*.id =*/"read3",
                        /*.sequence =*/"GGAGTATA"_dna5,
                        /*.base_qualities =*/"!!*+,-./"_phred42,
                        /*.alignment =*/alignment_t{"T-G-A-TC"_gapped_dna5, "G-AGTA-T"_gapped_dna5},
                        /*.reference_id =*/0, // "ref"
                        /*.reference_position =*/2,
                        /*.header_ptr =*/&header,
                        /*.flag =*/seqan3::sam_flag{43u},
                        {/*.mate_reference_id =*/0, /*.mate_position =*/9, /*.template_length =*/300},
                        /*.mapping_quality =*/63u,
                        /*.cigar_sequence =*/"1S1M1P1M1I1M1I1D1M1S"_cigar,
                        /*.tags =*/seqan3::sam_tag_dictionary{}};

    std::vector<record_type> records{record1, record2, record3};
};

} // namespace seqan3::test::fixture::io::sam_file

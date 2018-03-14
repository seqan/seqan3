// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ============================================================================

/*!\file
 * \brief Provides seqan3::sequence_file_in and corresponding traits classes.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <cassert>
#include <string>
#include <variant>
#include <vector>

#include <seqan3/alphabet/adaptation/all.hpp>
#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/range/concept.hpp>
#include <seqan3/range/container/concept.hpp>
#include <seqan3/range/container/concatenated_sequences.hpp>
#include <seqan3/io/sequence/sequence_file_format.hpp>
#include <seqan3/io/sequence/sequence_file_format_fasta.hpp>

namespace seqan3
{


template <typename t>
concept bool sequence_file_traits_concept = requires (t v)
{
    typename t::format_type;
    requires detail::all_holdees_meets_sequence_file_format_concept<typename t::format_type>();

    requires alphabet_concept<typename t::sequence_alphabet>;
    requires alphabet_concept<typename t::sequence_legal_alphabet>;
    requires sequence_concept<typename t::template sequence_container<typename t::sequence_alphabet>>;
    requires sequence_concept<typename t::template sequence_container_container<
        typename t::template sequence_container<typename t::sequence_alphabet>>>;

    requires alphabet_concept<typename t::id_alphabet>;
    requires sequence_concept<typename t::template id_container<typename t::id_alphabet>>;
    requires sequence_concept<typename t::template id_container_container<typename t::template id_container<
        typename t::id_alphabet>>>;

//     requires alphabet_concept<typename t::quality_alphabet>;
//     requires alphabet_concept<typename t::quality_legal_alphabet>;
//     requires sequence_concept<typename t::template quality_container<typename t::quality_alphabet>>;
//     requires sequence_concept<typename t::template quality_container_container<
//         typename t::template quality_container<typename t::quality_alphabet>>>;
};

struct sequence_file_default_traits_dna
{
    using format_type = std::variant<sequence_file_format_fasta>;
//                                      sequence_file_format_fastq<stream_type>,
//                                      sequence_file_format_embl<stream_type>,
//                                      sequence_file_format_genbank<stream_type>,
//                                      sequence_file_format_raw<stream_type>>;

    using sequence_alphabet                 = dna5;
    using sequence_legal_alphabet           = dna15;
    template <typename _sequence_alphabet>
    using sequence_container                = std::vector<_sequence_alphabet>;
    template <typename _sequence_container>
    using sequence_container_container      = concatenated_sequences<_sequence_container>;

    using id_alphabet                       = char;
    template <typename _id_alphabet>
    using id_container                      = std::basic_string<_id_alphabet>;
    template <typename _id_container>
    using id_container_container            = concatenated_sequences<_id_container>;

    using quality_alphabet                  = char; //TODO phred42;
    template <typename _quality_alphabet>
    using quality_container                 = std::vector<_quality_alphabet>;
    template <typename _quality_container>
    using quality_container_container       = concatenated_sequences<_quality_container>;
};

struct sequence_file_default_traits_aa : sequence_file_default_traits_dna
{
    using sequence_alphabet = aa27;
    using legal_alphabet = aa27;
};

} // namespace seqan3


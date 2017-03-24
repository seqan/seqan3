// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
//
// Copyright (c) 2006-2017, Knut Reinert, FU Berlin
// Copyright (c) 2016-2017, Knut Reinert & MPI Molekulare Genetik
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
// ==========================================================================
// Authors: Svenja Mehringer, Temesgen H. Dadi, Jongkyu Kim
//          <svenja.mehringer@fu-berlin.de> <temesgen.dadi@fu-berlin.de>
// ==========================================================================

#pragma once

#include <string>
#include <variant>
#include <vector>
#include <fstream>
#include <functional>

#include "../../container/concepts.hpp"
#include "../detail/file_base.hpp"
#include "sequence_file_format.hpp"
#include "sequence_file_format_fasta.hpp"

namespace seqan3
{

// ==================================================================
// sequence_file_in_traits
// ==================================================================

template <typename t>
concept bool sequence_file_in_traits_concept = requires (t v)
{
    typename t::stream_type;
    typename t::valid_format_types;
    requires detail::meets_sequence_file_format_concept<typename t::valid_format_types>(
        std::make_index_sequence<std::variant_size_v<typename t::valid_format_types>>{}
    );
    t::valid_compression_formats;
};

struct sequence_file_in_default_traits
{
    using stream_type = std::ifstream;
    using valid_format_types = std::variant<sequence_file_format_fasta/*,
                                       sequence_file_format_fastq,
                                       sequence_file_format_embl,
                                       sequence_file_format_genbank,
                                       sequence_file_format_raw*/>;
    using valid_compressions = std::variant<char/*specify compression formats*/>;
    static inline std::vector<std::pair<std::string, valid_compressions>> valid_compression_formats{};
};

// ==================================================================
// sequence_file_in
// ==================================================================

template <typename sequence_file_in_traits = sequence_file_in_default_traits>
    requires sequence_file_in_traits_concept<sequence_file_in_traits>
class sequence_file_in : public detail::file_base<sequence_file_in_traits>
{
public:
    /* types */
    using stream_type = typename detail::file_base<sequence_file_in_traits>::stream_type;
    using valid_format_types = typename detail::file_base<sequence_file_in_traits>::valid_format_types;

    /* constructors */
    sequence_file_in(std::experimental::filesystem::path _file_name) :
        detail::file_base<sequence_file_in_traits>(_file_name) {};
    sequence_file_in() = delete;
    sequence_file_in(sequence_file_in const &) = delete;
    sequence_file_in & operator=(sequence_file_in const &) = delete;
    sequence_file_in(sequence_file_in &&) = default;
    sequence_file_in & operator=(sequence_file_in &&) = default;
    ~sequence_file_in() = default;

    /* options */
    struct options_type
    {
        // post-processing filters that operate on buffer before assignment to out-value
        std::function<void(std::string &)> sequence_filter = [] (std::string & seq) {};
        std::function<void(std::string &)> meta_filter = [] (std::string & meta) {};
        std::function<void(std::string &)> qual_filter = [] (std::string & qual) {};
    };
    options_type options;

    /* member functions */
    // TODO make the requirements stricter
    template <typename sequence_type, typename meta_type, typename qual_type>
        requires sequence_concept<std::decay_t<sequence_type>> &&
                 sequence_concept<std::decay_t<meta_type>> &&
                 sequence_concept<std::decay_t<qual_type>>
    inline void read(sequence_type && seq,
                     meta_type && meta = std::string{},
                     qual_type && qual = std::string{});

    template <typename seqs_type, typename metas_type, typename quals_type>
        requires sequence_of_sequence_concept<std::decay_t<seqs_type>> &&
                 sequence_of_sequence_concept<std::decay_t<metas_type>> &&
                 sequence_of_sequence_concept<std::decay_t<quals_type>>
    inline void read(seqs_type && seqs,
                     metas_type && metas = std::vector<std::string>{},
                     quals_type && quals = std::vector<std::string>{},
                     size_t max_records = 0);

protected:
    /* member variables */
    using detail::file_base<sequence_file_in_traits>::format;
    using detail::file_base<sequence_file_in_traits>::stream;
};

// ------------------------------------------------------------------
// public API
// ------------------------------------------------------------------

template <typename sequence_file_in_traits>
template <typename sequence_type, typename meta_type, typename qual_type>
    requires sequence_concept<std::decay_t<sequence_type>> &&
             sequence_concept<std::decay_t<meta_type>> &&
             sequence_concept<std::decay_t<qual_type>>
inline void sequence_file_in<sequence_file_in_traits>::read(sequence_type && seq,
                                                            meta_type && meta,
                                                            qual_type && qual)
{
    if (format.valueless_by_exception())
        throw "format not set!";

    std::visit([&] (sequence_file_format_concept & f)
    {
        f.read(seq, meta, qual, stream, options);
    }, format);
}

template <typename sequence_file_in_traits>
template <typename seqs_type, typename metas_type, typename quals_type>
    requires sequence_of_sequence_concept<std::decay_t<seqs_type>> &&
             sequence_of_sequence_concept<std::decay_t<metas_type>> &&
             sequence_of_sequence_concept<std::decay_t<quals_type>>
inline void sequence_file_in<sequence_file_in_traits>::read(seqs_type && seqs,
                                                            metas_type && metas,
                                                            quals_type && quals,
                                                            size_t max_records)
{
    if (format.valueless_by_exception())
        throw "format not set!";

    std::visit([&] (sequence_file_format_concept & f)
    {
        f.read(seqs, metas, quals, stream, options, max_records);
    }, format);
}

} // namespace seqan

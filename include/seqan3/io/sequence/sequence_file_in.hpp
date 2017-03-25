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
// Author: Svenja Mehringer <svenja.mehringer@fu-berlin.de>
// Author: Temesgen H. Dadi <temesgen.dadi@fu-berlin.de>
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

namespace detail
{

// ==================================================================
// sequence_file_in_traits
// ==================================================================

template <typename type>
constexpr bool is_valid_compression_type(type const &)
{
    return false;
}

template <typename type>
constexpr bool is_valid_compression_type(std::vector<std::pair<std::string, type>> const &)
{
    return true;
}

//! A Concept that a sequence_file_in_traits object must satisfy
/*! When you want to instantiate a `sequence_file_in` object with
 * your own traits specification they must satisfy this concept.
 */
template <typename t>
concept bool sequence_file_in_traits_concept = requires (t v)
{
    typename t::stream_type;
    typename t::valid_format_types;

    requires detail::meets_sequence_file_format_concept<typename t::valid_format_types>(
        std::make_index_sequence<std::variant_size_v<typename t::valid_format_types>>{}
    );

    t::valid_compression_formats;
    requires is_valid_compression_type(t::valid_compression_formats);
};

} // namespace detail

// ==================================================================
// sequence_file_in_default_traits
// ==================================================================

//! The default configuration of seqan3::sequence_file_in
/*!
 * This implies the equivalence of
 * `seqan3::sequence_file_in{"myfile.fa"}` and
 * `seqan3::sequence_file_in<seqan3::sequence_file_in_default_traits>{"myfile.fa"}`.
 */
struct sequence_file_in_default_traits
{
    using stream_type = std::ifstream; //!< The default stream type is std::ifstream
    //! The default format types of the sequence file
    /*!
     * The valid types are stored as an `std::variant` of formats where each format
     * satisfies the sequence_file_format_concept.
     * Valid sequence formats by default are:
     * [fasta](https://www.genomatix.de/online_help/help/sequence_formats.html#FASTA),
     * [fastq](https://www.genomatix.de/online_help/help/sequence_formats.html#FASTQ),
     * [embl](https://www.genomatix.de/online_help/help/sequence_formats.html#EMBL),
     * [genebank](https://www.genomatix.de/online_help/help/sequence_formats.html#GB),
     * and
     * [raw](https://www.genomatix.de/online_help/help/sequence_formats.html#plain).
     */
    using valid_format_types = std::variant<sequence_file_format_fasta/*,
                                       sequence_file_format_fastq,
                                       sequence_file_format_embl,
                                       sequence_file_format_genbank,
                                       sequence_file_format_raw*/>;
    //! The default compression types that are supported
    /*!
      The compression types are stored as an `std::variant` of formats, which by
      default depend on the availability of external libraries like ZLIB.
    */
    using valid_compressions = std::variant<decltype(std::ignore)>;
    //! static member variable that stores valid compression formats
    /*!
      The valid_compression_formats stores pair of a compression format and its
      corresponding file extension as an identifier. Example: {".gz", compression_format_zlib}.
    */
    static inline std::vector<std::pair<std::string, valid_compressions>> valid_compression_formats{};
};

// ==================================================================
// sequence_file_in
// ==================================================================

//! A class that features reading sequence files from a stream
/*!
 * Use an instantiation of this class to read from a stream. You can specialize
 * the `sequence_file_in` object by specifying a `sequence_file_in_traits` as an
 * template argument otherwise the class is defaulted with `sequence_file_in_default_traits`.
 *
 * Example:
 *  \code{.cpp}
 *
 * int main()
 * {
 *     dna4_string seq;
 *     std::string id;
 *
 *     seqan3::sequence_file_in file_in{"input.fasta"};
 *
 *     file_in.read(seq, id);
 * }
 * \endcode
 *
 * \sa sequence_file_in_default_traits
 */
template <typename sequence_file_in_traits = sequence_file_in_default_traits>
    requires detail::sequence_file_in_traits_concept<sequence_file_in_traits>
class sequence_file_in : protected detail::file_base<sequence_file_in_traits>
{
public:
    using detail::file_base<sequence_file_in_traits>::stream_type;
    using detail::file_base<sequence_file_in_traits>::valid_format_types;
    /* constructors */
    //! constructor with file name argument
    /*!
     * Passing a file name (path) as an argument to the constructor will open
     * the stream using this name.
     *
     * Note: The sequence file format will automatically be deduced by
     * the extension file name extension:
     *
     *    -# Check whether a valid compression format was used. Is so, strip the
     *        strip the file name and continue, if not continue with the original
     *        file name
     *
     *    -# Check every valid file format in `valid_format_types` for their
     *        extension identifiers in `file_extensions` and choose the according
     *        format. this function will __throw__ when the format cannot be inferred.
     */
    explicit sequence_file_in(std::experimental::filesystem::path _file_name) :
        detail::file_base<sequence_file_in_traits>(std::move(_file_name)) {};
    sequence_file_in() = delete;
    sequence_file_in(sequence_file_in const &) = delete;
    sequence_file_in & operator=(sequence_file_in const &) = delete;
    sequence_file_in(sequence_file_in &&) = default; //!< default move constructor
    sequence_file_in & operator=(sequence_file_in &&) = default; //!< default move assignment constructor
    ~sequence_file_in() = default; //!< default deconstructor

    //! A struct holding additional features for the sequence_file_in object
    /*!
     * The options_type struct stores three std::functions that can alter the
     * sequence information directly while reading. The default functions do not
     * change the input but you can assign different function that crop, replace
     * or append to the sequence information (sequence_filter), the meta
     * information (meta_filter) or the quality information (qual_filter).
     *
     * For example this code snippet will only extract the first 5 characters
     * of the sequence identifier (meta information):
     *
     * Example:
     *
     * \code{.cpp}
     * int main()
     * {
     *      dna4_string seq;
     *      std::string id;
     *
     *      seqan3::sequence_file_in file_in{"input.fasta"};
     *
     *      // crop the meta information
     *      file_in.options.meta_filter = [](std::string & in){return in.substr(0, 5);};
     *
     *      file_in.read(seq, id);
     * }
     * \endcode
     */
    struct options_type
    {
        // post-processing filters that operate on buffer before assignment to out-value
        std::function<void(std::string &)> sequence_filter = [] (std::string & seq) {}; //!< alters the raw sequence
        std::function<void(std::string &)> meta_filter = [] (std::string & meta) {}; //!< alters meta information
        std::function<void(std::string &)> qual_filter = [] (std::string & qual) {}; //!< alters the quality sequence
    };
    options_type options; //!< holds the filter functions

    // TODO make the requirements stricter
    //! reads a single record from the stream and into the given arguments
    /*!
     * \param seq the raw sequence information.
     * \param meta the meta information (e.g. the sequence identifier/name).
     * \param qual the quality information.
     */
    template <typename sequence_type, typename meta_type, typename qual_type>
        requires sequence_concept<std::decay_t<sequence_type>> &&
                 sequence_concept<std::decay_t<meta_type>> &&
                 sequence_concept<std::decay_t<qual_type>>
    void read(sequence_type && seq,
              meta_type && meta = std::string{},
              qual_type && qual = std::string{});

    // TODO make the requirements stricter
    //! reads many or all information from the stream appending it to the given arguments
    /*!
     * \param seqs a container of sequences to append to.
     * \param metas a container of meta information to append to.
     * \param quals a container of quality information.
     * \param max_records limit the number of records to read to max_records.
     */
    template <typename seqs_type, typename metas_type, typename quals_type>
        requires sequence_of_sequence_concept<std::decay_t<seqs_type>> &&
                 sequence_of_sequence_concept<std::decay_t<metas_type>> &&
                 sequence_of_sequence_concept<std::decay_t<quals_type>>
    void read(seqs_type && seqs,
              metas_type && metas = std::vector<std::string>{},
              quals_type && quals = std::vector<std::string>{},
              size_t max_records = 0);

protected:
    /* member variables */
    using detail::file_base<sequence_file_in_traits>::format;
    using detail::file_base<sequence_file_in_traits>::stream;
};

// ------------------------------------------------------------------
// public functions
// ------------------------------------------------------------------

template <typename sequence_file_in_traits>
template <typename sequence_type, typename meta_type, typename qual_type>
    requires sequence_concept<std::decay_t<sequence_type>> &&
             sequence_concept<std::decay_t<meta_type>> &&
             sequence_concept<std::decay_t<qual_type>>
inline void
sequence_file_in<sequence_file_in_traits>::read(sequence_type && seq,
                                                meta_type && meta,
                                                qual_type && qual)
{
    if (format.valueless_by_exception())
        throw std::runtime_error("format not set!"); // TODO:: replace by internal seqan3 error

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
inline void
sequence_file_in<sequence_file_in_traits>::read(seqs_type && seqs,
                                                metas_type && metas,
                                                quals_type && quals,
                                                size_t max_records)
{
    if (format.valueless_by_exception())
        throw std::runtime_error("format not set!"); // TODO:: replace by internal seqan3 error

    std::visit([&] (sequence_file_format_concept & f)
    {
        f.read(seqs, metas, quals, stream, options, max_records);
    }, format);
}

} // namespace seqan

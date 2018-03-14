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
#include <fstream>
#include <string>
#include <variant>
#include <vector>

#include <seqan3/io/record.hpp>
#include <seqan3/io/sequence/sequence_file_traits.hpp>
#include <seqan3/io/detail/in_file_iterator.hpp>
#include <seqan3/io/sequence/sequence_file_format.hpp>
#include <seqan3/io/sequence/sequence_file_format_fasta.hpp>


namespace seqan3::detail
{

// template <typename type>
// concept bool sequence_file_in_valid_fields = requires (type const & val)
// {
//     for (field f : type::as_array)
//         require (f == field::SEQ || f == field::ID || f == field::QUAL || f == field::SEQ_QUAL);
// };

// template <typename tup_t = std::tuple<>, size_t field_no = 0, typename ...types, field ...fs>
// constexpr auto generate_record_type(meta::list<types...> const &, fields<fs...> const &)
// {
//     if constexpr (field_no == ffields<fs...>::as_array.size())
//         return tup_t{};
//
//     using new_tup_t = decltype(std::tuple_cat(tup{},
//                                               std::tuple<meta::at<meta::list<types...>,
//                                                                   static_cast<size_t>(fields<fs...>::as_array[field_no])>{}));
//     return generate_record_type<new_tup_t, field_no++>(
// }
//


}

namespace seqan3
{

/*!\brief A class for reading sequence files, e.g. FASTA, FASTQ ...
 * \ingroup io
 * \tparam traits_type An auxiliary type that defines certain member type and constants, must satisfy
 * seqan3::sequence_file_traits_concept.
 * \tparam stream_type The type of the stream, must satisfy seqan3::istream_concept.
 * \tparam selected_fields A \ref fields type with the list and order of desired record entries.
 *
 * \details
 *
 * ### Introduction general
 *
 * TODO move this to io/all.hpp
 *
 * ### Introduction
 *
 * Sequence files are the most generic and common biological files. Well-known formats include
 * FastA and FastQ, but some may also be interested in treating SAM or BAM files as sequence
 * files, discarding the alignment.
 *
 * The Sequence file abstraction provides two fields: seqan3::field::SEQ and seqan3::field::ID.
 *
 * By default the SEQ field is retrieved as a vector over seqan3::qualified <seqan3::dna5>, i.e.
 * the sequence combines qualities and actual sequence in one. You can later drop the qualities if
 * you are not interested in them, or specify a custom traits type to change the underlying
 * alphabet of the sequence field so they are never returned.
 *
 * ### Construction and specialisation
 *
 * This class comes with two constructors, one for construction from a file name and one for construction from
 * an existing stream and a known format. The first one automatically picks the format based on the extension
 * of the file name. The second can be used if you have a non-file stream, like std::cin or std::istringstream,
 * that you want to read from and/or if you cannot use file-extension based detection, but know that your input
 * file has a certain format.
 *
 * In most cases the template parameters are deduced completely automatically:
 *
 * ```cpp
 * sequence_file_in fin{"/tmp/my.fasta"}; // FastA with DNA sequences assumed, regular std::ifstream taken as stream
 * ```
 * Reading from an std::istringstream:
 * ```cpp
 * std::string input
 * {
 *     "> TEST1\n"
 *     "ACGT\n"
 *     "> Test2\n"
 *     "AGGCTGN\n"
 *     "> Test3\n"
 *     "GGAGTATAATATATATATATATAT\n"
 * };
 *
 * std::istringstream iss(input);
 *
 * sequence_file_in fin{std::move(iss), sequence_file_format_fasta{}};
 * //              ^ no need to specify the template arguments
 * ```
 *
 * In some cases, you do need to specify the arguments, e.g. if you want to read amino acids:
 *
 * ```cpp
 * sequence_file_in<seqeuence_file_default_traits_aa> fin{"/tmp/my.fasta"};
 * ```
 *
 * You can define your own traits type to further customise the types used by and return by this class, see
 * seqan3::seqeuence_file_default_traits_dna for more details.
 *
 * ### Reading record-wise
 *
 * You can iterate over this file record-wise:
 *
 * ```cpp
 * sequence_file_in fin{"/tmp/my.fasta"};
 *
 * for (auto record : fin)
 * {
 *     std::cout << "ID:  " << std::get<1>(record) << '\n';
 *     std::cout << "SEQ: " << (std::get<0>(record) | view::to_char) << '\n'; // sequence is converted to char string on-the-fly
 * }
 * ```
 *
 * In the above example, record has the type \ref record_type which is a tuple, that's why we can access it via
 * std::get.
 *
 * However we can also directly use [structured bindings](http://en.cppreference.com/w/cpp/language/structured_binding)
 * to decompose the tuple into its elements:
 *
 * ```cpp
 * sequence_file_in fin{"/tmp/my.fasta"};
 *
 * for (auto [ seq, id ] : fin)
 * {
 *     std::cout << "ID:  " << id << '\n';
 *     std::cout << "SEQ: " << (seq | view::to_char) << '\n'; // sequence string is converted to char string on-the-fly
 * }
 * ```
 * In this case you immediately get the two elements of the tuple: `seq` of \ref sequence_type and `id` of
 * \ref id_type.
 *
 * ### Reading column-wise
 *
 * TODO example
 *
 *
 * ### Formats
 *
 * TODO give overview of formats
 *
 *
 *
 */
template <sequence_file_traits_concept traits_type = sequence_file_default_traits_dna,
          typename stream_type = std::ifstream,
          typename selected_fields = fields<field::SEQ, field::ID>> // TODO istream_concept stream_type
class sequence_file_in
{
public:
    using format_type = typename traits_type::format_type;

    /* record and column types */
    /*!\name Data member types
     * \{
     */
    //!\brief The type of the sequence field (usually std::vector of seqan3::dna5q).
    using sequence_type         = typename traits_type::template sequence_container<typename traits_type::sequence_alphabet>;
    //!\brief The type of the ID field (usually std::string).
    using id_type               = typename traits_type::template id_container<typename traits_type::id_alphabet>;
    //!\brief The type of the QUAL field TODO
    using quality_type             = typename traits_type::template quality_container<typename traits_type::quality_alphabet>;
////!\brief The type of the QUAL field TODO
//     using seq_quality_type         = typename traits_type::
//                                     template sequence_container<quality_composition<typename traits_type::sequence_alphabet,
//                                                                                     typename traits_type::quality_alphabet>>;
    //!\brief The previously defined types in a type list (meta::list).
    using field_types           = meta::list<sequence_type, id_type, quality_type, quality_type>; //TODO replace last with seq_quality_type
    //!\brief The corresponding seqan3::field IDs to \ref field_types.
    using field_types_as_ids    = fields<field::SEQ, field::ID, field::QUAL, field::SEQ_QUAL>;
    //!\brief The type of the record, specialisation of seqan3::record; acts as a tuple of the selected field types.
    using record_type           = record<field_types, field_types_as_ids, selected_fields>;
    //!\brief Column type of field 0 (sequences); usually seqan3::concatenated_sequences <sequence_type>.
    using sequence_column_type  = typename traits_type::template sequence_container_container<sequence_type>;
    //!\brief Column type of field 1 (IDs); usually seqan3::concatenated_sequences <id_type>.
    using id_column_type        = typename traits_type::template id_container_container<id_type>;
    //!\brief A tuple of sean3::field that specifies the order of the fields.
//     static constexpr std::tuple<field, field> fields = { field::SEQ, field::ID };
    //!\}

    /*!\name Range member types
     * \{
     */
    //!\brief The value_type is the recored_type.
    using value_type        = record_type;
    //!\brief The reference_type.
    using reference         = record_type &;
    //!\brief The const_reference type.
    using const_reference   = record_type const &;
    //!\brief An unsigned integer type, usually std::size_t.
    using size_type         = size_t;
    //!\brief A signed integer type, usually std::ptrdiff_t.
    using difference_type   = std::make_signed_t<size_t>;
    //!\brief The iterator type of this view (an input iterator).
    using iterator          = detail::in_file_iterator<sequence_file_in>;
    //!\brief The const iterator has an invalid type. TODO document that is not const-iterable
    using const_iterator    = detail::in_file_iterator<sequence_file_in const>;
    //!\brief The iterator type of this view (an input iterator).
    using sentinel          = detail::in_file_sentinel<sequence_file_in>;
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    //!\brief Construct from filename.
    sequence_file_in(std::string const & _file_name); //TODO fs::path

    /*!\brief Construct from an existing stream and with specified format.
     * \tparam file_format The format of the file in the stream, must satisfy seqan3::sequence_file_format_concept.
     * \param[in] _stream The stream to operate on (this must be std::move'd in!).
     * \param[in] tag The file format tag.
     *
     * \details
     *
     * TODO
     */
    template <sequence_file_format_concept file_format>
    sequence_file_in(stream_type && _stream, file_format const & tag);

    //!\brief Copy construction is explicitly deleted, because you can't have multiple access to the same file.
    sequence_file_in(sequence_file_in const &) = delete;
    //!\brief Copy assignment is explicitly deleted, because you can't have multiple access to the same file.
    sequence_file_in & operator=(sequence_file_in const &) = delete;

    //!\brief Move construction is defaulted.
    sequence_file_in(sequence_file_in &&) = default;
    //!\brief Move assignment is defaulted.
    sequence_file_in & operator=(sequence_file_in &&) = default;
    //!\brief Destructor is defaulted.
    ~sequence_file_in() = default;
    //!\}

    /*!\name Range interface
     * \brief Provides functions for record based reading of the file.
     * \{
     */
    /*!\brief Returns an iterator to current position in the file.
     * \returns An iterator pointing to the current position in the file.
     *
     * Equals end() if the file is at end.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    iterator begin() noexcept
    {
        return {*this};
    }

    //!\brief TODO document impotence
    const_iterator cbegin() const noexcept;

    /*!\brief Returns a sentinel for comparison with iterator.
     * \returns Iterator to the first element.
     *
     * This element acts as a placeholder; attempting to dereference it results in undefined behaviour.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    sentinel end() noexcept
    {
        return {};
    }

    //!\copydoc end
    sentinel cend() const noexcept;

    /*!\brief Return the record we are currently at in the file.
     * \returns The current record marked as an rvalue.
     *
     * The record is moved out of the file upon calling this so you avoid expensive copies. This
     * also means that if you call this function a second time without incrementing an iterator,
     * you will receive an empty record.
     *
     * ### Example
     *
     * ```cpp
     *
     * sequence_file_in fin{"/tmp/my.fasta"};
     *
     * auto record1 = fin.back();         // record is retrieved and moved into record1
     * auto record2 = fin.back();         // the buffer is now empty so record2 received nothing
     *
     * auto it = fin.begin();
     *
     * ++it;                              // file iterator is moved to next record
     * auto record3 = fin.back();         // record3 now receives data again
     *
     * ++it;                              // file iterator is moved to next record
     * auto const & record4 = fin.back(); // this binds to the reference and does not trigger a move
     * auto const & record5 = fin.back(); // of course you can have multiple references
     * ++it;                              // however the previous references may now be invalid
     * ```
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    reference back() noexcept
    {
        return record_buffer;
    }

    //!\copydoc back
    const_reference back() const noexcept
    {
        return record_buffer;
    }
    //!\}


    /*!\name Tuple interface
     * \brief Provides functions for field-based ("column"-based) reading.
     * \{
     */
    template <size_t i>
    friend auto & get(sequence_file_in & file)
    {
        static_assert(i < 2);

        file.read_columns();

        if constexpr (i == 0)
            return file.sequence_column_buffer;
        else
            return file.id_column_buffer;
    }

    template <size_t i>
    friend auto && get(sequence_file_in && file)
    {
        return std::move(get<i>(file));
    }

//     template <field f>
//     friend auto && get(sequence_file_in & file)
//     {
//         return get<get_position_in_tuple(f, fields)>(file);
//     }
    //!\}

    /* options */
    struct options_type
    {
        // post-processing filters that operate on buffer before assignment to out-value
        std::function<void(std::string &)>  sequence_filter = [] (std::string &) {};
        std::function<void(std::string &)>        id_filter = [] (std::string &) {};
        std::function<void(std::string &)>      quality_filter = [] (std::string &) {};
    };
    options_type options;

private:
    /* buffers */
    record_type             record_buffer;
    sequence_column_type    sequence_column_buffer;
    id_column_type          id_column_buffer;

    /* file format */
    std::string file_name;
    stream_type stream;
    format_type format;

    /* private functions */
//     void select_decompression(std::string const & compress_ext);
    template <size_t index = 0>
    void select_format(std::string const & ext);

    void read_record(record_type &);
    void buffer_next_record();
    void read_columns();

    // befriend iterator so it can access record, stream, and read()...
    friend iterator;
    friend const_iterator;
};

template <sequence_file_traits_concept traits_type, typename stream_type, typename selected_fields>
sequence_file_in(std::string const &)
    -> sequence_file_in<sequence_file_default_traits_dna,
                        std::ifstream,
                        fields<field::SEQ, field::ID>>;// TODO virtual_stream<input>

template <sequence_file_traits_concept traits_type, typename stream_type, typename selected_fields, sequence_file_format_concept file_format>
sequence_file_in(stream_type && _stream, file_format const &)
    -> sequence_file_in<sequence_file_default_traits_dna,
                        std::decay_t<stream_type>,
                        fields<field::SEQ, field::ID>>;


// ------------------------------------------------------------------
// public API
// ------------------------------------------------------------------

template <sequence_file_traits_concept traits_type, typename stream_type, typename selected_fields>
sequence_file_in<traits_type, stream_type, selected_fields>::sequence_file_in(std::string const & _file_name) :
        file_name(_file_name)
{
//     // open stream
//     stream.open(_file_name, std::ios::binary);
//
//     // initialise format handler
//     std::string ext{get_file_extension(file_name)};
//     select_format<0>(format, ext);
}

template <sequence_file_traits_concept traits_type, typename stream_type, typename selected_fields>
template <sequence_file_format_concept file_format>
sequence_file_in<traits_type, stream_type, selected_fields>::sequence_file_in(stream_type && _stream, file_format const &) :
    stream{std::move(_stream)}, format{file_format{}}
{
    buffer_next_record();
}

template <sequence_file_traits_concept traits_type, typename stream_type, typename selected_fields>
inline void
sequence_file_in<traits_type, stream_type, selected_fields>::read_record(record_type & _record)
{
    assert(!format.valueless_by_exception);
    std::visit([&] (sequence_file_format_concept & f) { f.read(_record, stream, options); }, format);
}

template <sequence_file_traits_concept traits_type, typename stream_type, typename selected_fields>
inline void
sequence_file_in<traits_type, stream_type, selected_fields>::buffer_next_record()
{
    read_record(record_buffer);
}

template <sequence_file_traits_concept traits_type, typename stream_type, typename selected_fields>
inline void
sequence_file_in<traits_type, stream_type, selected_fields>::read_columns()
{
    //TODO delegate to format if column-wise available there
    if (!stream.eof())
    {
        // read the remaining records and split into column buffers
        for (auto [ seq, id ] : std::move(*this))
        {
            sequence_column_buffer.push_back(std::move(seq));
            id_column_buffer.push_back(std::move(id));
        }
    }
}

// ------------------------------------------------------------------
// private functions
// ------------------------------------------------------------------

// inline void
// sequence_file_in::select_decompression(std::string const & compress_ext)
// {
//     for (auto const & pair : valid_compression_formats)
//     {
//         if (compress_ext == std::get<0>(pair))
//         {
//             std::visit([&stream] (auto & compressor) { stream.push(compressor); }, std::get<1>(pair));
//             break;
//         }
//     }
// }

template <sequence_file_traits_concept traits_type, typename stream_type, typename selected_fields>
template <size_t index>
inline void
sequence_file_in<traits_type, stream_type, selected_fields>::select_format(std::string const & ext)
{
    if (index == std::variant_size_v<format_type>)
        throw std::runtime_error("No valid format found for this extension.");
    else if (std::variant_alternative_t<index, format_type>::file_extensions().contains(ext))
        format = std::variant_alternative_t<index, format_type>{};
    else
        select_format<index+1>(format, ext);
}

} // namespace seqan3


// ------------------------------------------------------------------
// public API
// ------------------------------------------------------------------


namespace std
{

template <typename traits_type, typename stream_type, typename selected_fields>
struct tuple_size<seqan3::sequence_file_in<traits_type, stream_type, selected_fields>>
{
    static constexpr size_t value = 2;
};

template <typename traits_type, typename stream_type, typename selected_fields>
struct tuple_element<0, seqan3::sequence_file_in<traits_type, stream_type, selected_fields>>
{
    using type = typename seqan3::sequence_file_in<traits_type, stream_type, selected_fields>::sequence_column_type;
};

template <typename traits_type, typename stream_type, typename selected_fields>
struct tuple_element<1, seqan3::sequence_file_in<traits_type, stream_type, selected_fields>>
{
    using type = typename seqan3::sequence_file_in<traits_type, stream_type, selected_fields>::id_column_type;
};

}

/*
template <seqan3::field f, typename ...types>
auto get(sequence_file_in<types...>::record_type & tup)
{
    static constexpr size_t i = [] ()
    {
        for (size_t j = 0; j < std::tuple_size_v<sequence_file_in<types...>::record_type>; ++j)
            if (


    };



}
*/

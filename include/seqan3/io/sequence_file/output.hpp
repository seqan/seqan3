// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::sequence_file_output and corresponding traits classes.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <cassert>
#include <filesystem>
#include <fstream>
#include <ranges>
#include <string>
#include <variant>
#include <vector>

#include <seqan3/io/detail/misc.hpp>
#include <seqan3/io/detail/misc_output.hpp>
#include <seqan3/io/detail/out_file_iterator.hpp>
#include <seqan3/io/detail/record.hpp>
#include <seqan3/io/detail/record_like.hpp>
#include <seqan3/io/exception.hpp>
#include <seqan3/io/record.hpp>
#include <seqan3/io/sam_file/format_sam.hpp>
#include <seqan3/io/sequence_file/format_embl.hpp>
#include <seqan3/io/sequence_file/format_fasta.hpp>
#include <seqan3/io/sequence_file/format_fastq.hpp>
#include <seqan3/io/sequence_file/format_genbank.hpp>
#include <seqan3/io/sequence_file/output_format_concept.hpp>
#include <seqan3/io/sequence_file/output_options.hpp>
#include <seqan3/io/stream/concept.hpp>
#include <seqan3/utility/tuple/concept.hpp>
#include <seqan3/utility/type_list/traits.hpp>
#include <seqan3/utility/views/convert.hpp>
#include <seqan3/utility/views/elements.hpp>
#include <seqan3/utility/views/zip.hpp>

namespace seqan3
{

// ----------------------------------------------------------------------------
// sequence_file_output
// ----------------------------------------------------------------------------

/*!\brief A class for writing sequence files, e.g. FASTA, FASTQ ...
 * \ingroup io_sequence_file
 * \tparam selected_field_ids   A seqan3::fields type with the list and order of fields IDs; only relevant if these
 *                              can't be deduced.
 * \tparam valid_formats        A seqan3::type_list of the selectable formats (each must meet
 *                              seqan3::sequence_file_output_format).
 *
 * \details
 *
 * \include{doc} doc/fragments/sequence_file_output.md
 *
 * \remark For a complete overview, take a look at \ref io_sequence_file
 */

template <detail::fields_specialisation selected_field_ids_ = fields<field::seq, field::id, field::qual>,
          detail::type_list_of_sequence_file_output_formats valid_formats_ =
              type_list<format_embl, format_fasta, format_fastq, format_genbank, format_sam>>
class sequence_file_output
{
public:
    /*!\name Template arguments
     * \brief Exposed as member types for public access.
     * \{
     */
    //!\brief A seqan3::fields list with the fields selected for the record.
    using selected_field_ids = selected_field_ids_;
    //!\brief A seqan3::type_list with the possible formats.
    using valid_formats = valid_formats_;
    //!\brief Character type of the stream(s).
    using stream_char_type = char;
    //!\}

    //!\brief The subset of seqan3::field IDs that are valid for this file.
    using field_ids = fields<field::seq, field::id, field::qual>;

    static_assert(
        []() constexpr
        {
            for (field f : selected_field_ids::as_array)
                if (!field_ids::contains(f))
                    return false;
            return true;
        }(),
        "You selected a field that is not valid for sequence files, please refer to the documentation "
        "of sequence_file_output::field_ids for the accepted values.");

    /*!\name Range associated types
     * \brief Most of the range associated types are `void` for output ranges.
     * \{
     */

    //!\brief The value type (void).
    using value_type = void;
    //!\brief The reference type (void).
    using reference = void;
    //!\brief The const reference type (void).
    using const_reference = void;
    //!\brief The size type (void).
    using size_type = void;
    //!\brief A signed integer type, usually std::ptrdiff_t.
    using difference_type = std::ptrdiff_t;
    //!\brief The iterator type of this view (an output iterator).
    using iterator = detail::out_file_iterator<sequence_file_output>;
    //!\brief The const iterator type is void, because files are not const-iterable.
    using const_iterator = void;
    //!\brief The type returned by end().
    using sentinel = std::default_sentinel_t;
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    //!\brief Default constructor is explicitly deleted, you need to give a stream or file name.
    sequence_file_output() = delete;
    //!\brief Copy construction is explicitly deleted, because you can't have multiple access to the same file.
    sequence_file_output(sequence_file_output const &) = delete;
    //!\brief Copy assignment is explicitly deleted, because you can't have multiple access to the same file.
    sequence_file_output & operator=(sequence_file_output const &) = delete;
    //!\brief Move construction is defaulted.
    sequence_file_output(sequence_file_output &&) = default;
    //!\brief Move assignment is defaulted.
    sequence_file_output & operator=(sequence_file_output &&) = default;
    //!\brief Destructor is defaulted.
    ~sequence_file_output() = default;

    /*!\brief Construct from filename.
     * \param[in] filename      Path to the file you wish to open.
     * \param[in] fields_tag    A seqan3::fields tag. [optional]
     *
     * \details
     *
     * In addition to the file name, you may specify a custom seqan3::fields type which may be easier than
     * defining all the template parameters.
     *
     * ### Compression
     *
     * This constructor transparently applies a compression stream on top of the file stream in case
     * the given file extension suggests the user wants this.
     * See the section on \link io_compression compression and decompression \endlink for more information.
     */
    sequence_file_output(std::filesystem::path filename,
                         selected_field_ids const & SEQAN3_DOXYGEN_ONLY(fields_tag) = selected_field_ids{}) :
        primary_stream{new std::ofstream{}, stream_deleter_default}
    {
        primary_stream->rdbuf()->pubsetbuf(stream_buffer.data(), stream_buffer.size());
        static_cast<std::basic_ofstream<char> *>(primary_stream.get())
            ->open(filename, std::ios_base::out | std::ios::binary);

        if (!primary_stream->good())
            throw file_open_error{"Could not open file " + filename.string() + " for writing."};

        // possibly add intermediate compression stream
        secondary_stream = detail::make_secondary_ostream(*primary_stream, filename);

        // initialise format handler or throw if format is not found
        detail::set_format(format, filename);
    }

    /*!\brief Construct from an existing stream and with specified format.
     * \tparam file_format   The format of the file in the stream, must satisfy seqan3::sequence_file_output_format.
     * \param[in,out] stream The stream to write to, must be derived of std::basic_ostream<stream_char_t>.
     * \param[in] format_tag The file format tag.
     * \param[in] fields_tag A seqan3::fields tag. [optional]
     *
     * \details
     *
     * ### Compression
     *
     * This constructor **does not** apply compression transparently (because there is no way to know if the user
     * wants this). However, you can just pass e.g. seqan3::contrib::gz_ostream to this constructor if you explicitly
     * want compression.
     * See the section on \link io_compression compression and decompression \endlink for more information.
     */
    template <output_stream stream_t, sequence_file_output_format file_format>
        requires std::same_as<typename std::remove_reference_t<stream_t>::char_type, stream_char_type>
    sequence_file_output(stream_t & stream,
                         file_format const & SEQAN3_DOXYGEN_ONLY(format_tag),
                         selected_field_ids const & SEQAN3_DOXYGEN_ONLY(fields_tag) = selected_field_ids{}) :
        primary_stream{&stream, stream_deleter_noop},
        secondary_stream{&stream, stream_deleter_noop},
        format{detail::sequence_file_output_format_exposer<file_format>{}}
    {
        static_assert(list_traits::contains<file_format, valid_formats>,
                      "You selected a format that is not in the valid_formats of this file.");
    }

    //!\overload
    template <output_stream stream_t, sequence_file_output_format file_format>
        requires std::same_as<typename std::remove_reference_t<stream_t>::char_type, stream_char_type>
    sequence_file_output(stream_t && stream,
                         file_format const & SEQAN3_DOXYGEN_ONLY(format_tag),
                         selected_field_ids const & SEQAN3_DOXYGEN_ONLY(fields_tag) = selected_field_ids{}) :
        primary_stream{new stream_t{std::move(stream)}, stream_deleter_default},
        secondary_stream{&*primary_stream, stream_deleter_noop},
        format{detail::sequence_file_output_format_exposer<file_format>{}}
    {
        static_assert(list_traits::contains<file_format, valid_formats>,
                      "You selected a format that is not in the valid_formats of this file.");
    }
    //!\}

    /*!\name Range interface
     * \brief Provides functions for record based writing of the file.
     * \{
     */
    /*!\brief Returns an iterator to current position in the file.
     * \returns An iterator pointing to the current position in the file.
     *
     * You can write to the file by assigning to the iterator, but using push_back() is usually more intuitive.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     *
     * ### Example
     *
     * \include test/snippet/io/sequence_file/sequence_file_output_range_interface.cpp
     */
    iterator begin() noexcept
    {
        return {*this};
    }

    /*!\brief Returns a sentinel for comparison with iterator.
     * \returns An end that is never reached.
     *
     * This element acts as a placeholder; attempting to dereference it results in undefined behaviour. It
     * always compares false against an iterator.
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

    /*!\brief           Write a seqan3::record to the file.
     * \tparam record_t Type of the record, a specialisation of seqan3::record.
     * \param[in] r     The record to write.
     *
     * \details
     *
     * ### Complexity
     *
     * Constant. TODO linear in the size of the written sequences?
     *
     * ### Exceptions
     *
     * Basic exception safety.
     *
     * ### Example
     *
     * \include test/snippet/io/sequence_file/sequence_file_output_push_back_record.cpp
     */
    template <typename record_t>
    void push_back(record_t && r)
        requires detail::record_like<record_t>
    {
        write_record(detail::get_or_ignore<field::seq>(r),
                     detail::get_or_ignore<field::id>(r),
                     detail::get_or_ignore<field::qual>(r));
    }

    /*!\brief           Write a record in form of a std::tuple to the file.
     * \tparam tuple_t  Type of the record, a specialisation of std::tuple.
     * \param[in] t     The record to write.
     *
     * \details
     *
     * The fields in the tuple are assumed to correspond to the field IDs given in selected_field_ids, however
     * passing less is accepted if the format does not require all of them.
     *
     * ### Complexity
     *
     * Constant. TODO linear in the size of the written sequences?
     *
     * ### Exceptions
     *
     * Basic exception safety.
     *
     * ### Example
     *
     * \include test/snippet/io/sequence_file/sequence_file_output_push_back_tuple.cpp
     */
    template <typename tuple_t>
    void push_back(tuple_t && t)
        requires tuple_like<tuple_t> && (!detail::record_like<tuple_t>)
    {
        // index_of might return npos, but this will be handled well by get_or_ignore (and just return ignore)
        write_record(detail::get_or_ignore<selected_field_ids::index_of(field::seq)>(t),
                     detail::get_or_ignore<selected_field_ids::index_of(field::id)>(t),
                     detail::get_or_ignore<selected_field_ids::index_of(field::qual)>(t));
    }

    /*!\brief            Write a record to the file by passing individual fields.
     * \tparam arg_t     Type of the first field.
     * \tparam arg_types Types of further fields.
     * \param[in] arg    The first field to write.
     * \param[in] args   Further fields.
     *
     * \details
     *
     * The fields are assumed to correspond to the field IDs given in selected_field_ids, however passing less
     * is accepted if the format does not require all of them.
     *
     * ### Complexity
     *
     * Constant. TODO linear in the size of the written sequences?
     *
     * ### Exceptions
     *
     * Basic exception safety.
     *
     * ### Example
     *
     * \include test/snippet/io/sequence_file/sequence_file_output_emplace_back.cpp
     */
    template <typename arg_t, typename... arg_types>
    void emplace_back(arg_t && arg, arg_types &&... args)
    {
        push_back(std::tie(arg, args...));
    }

    /*!\brief            Write a range of records (or tuples) to the file.
     * \tparam rng_t     Type of the range, must satisfy std::ranges::output_range and have a reference type that
     *                   satisfies seqan3::tuple_like.
     * \param[in] range  The range to write.
     *
     * \details
     *
     * This function simply iterates over the argument and calls push_back() on each element.
     *
     * ### Complexity
     *
     * Linear in the number of records.
     *
     * ### Exceptions
     *
     * Basic exception safety.
     *
     * ### Example
     *
     * \include test/snippet/io/sequence_file/sequence_file_output_batch_write.cpp
     */
    template <std::ranges::input_range rng_t>
    sequence_file_output & operator=(rng_t && range)
        requires tuple_like<std::ranges::range_reference_t<rng_t>>
    {
        for (auto && record : range)
            push_back(std::forward<decltype(record)>(record));
        return *this;
    }

    /*!\brief            Write a range of records (or tuples) to the file.
     * \tparam rng_t     Type of the range, must satisfy std::ranges::input_range and have a reference type that
     *                   satisfies seqan3::tuple_like.
     * \param[in] range  The range to write.
     * \param[in] f      The file being written to.
     *
     * \details
     *
     * This operator enables sequence_file_output to be at the end of a piping operation. It just calls
     * operator=() internally.
     *
     * ### Complexity
     *
     * Linear in the number of records.
     *
     * ### Exceptions
     *
     * Basic exception safety.
     *
     * ### Example
     *
     * \include test/snippet/io/sequence_file/sequence_file_output_batch_write.cpp
     *
     * This is especially useful in combination with file-based filters:
     *
     * \snippet test/snippet/io/sequence_file/sequence_file_output_view_pipeline.cpp snippet
     */
    template <std::ranges::input_range rng_t>
    friend sequence_file_output & operator|(rng_t && range, sequence_file_output & f)
        requires tuple_like<std::ranges::range_reference_t<rng_t>>
    {
        f = range;
        return f;
    }

    //!\overload
    template <std::ranges::input_range rng_t>
    friend sequence_file_output operator|(rng_t && range, sequence_file_output && f)
        requires tuple_like<std::ranges::range_reference_t<rng_t>>
    {
        f = range;
        return std::move(f);
    }
    //!\}

    //!\brief The options are public and its members can be set directly.
    sequence_file_output_options options{};

    /*!\cond DEV
     * \brief Expose a reference to the secondary stream. [public, but not documented as part of the API]
     */
    std::basic_ostream<stream_char_type> & get_stream()
    {
        return *secondary_stream;
    }
    //!\endcond

protected:
    //!\privatesection
    //!\brief A larger (compared to stl default) stream buffer to use when reading from a file.
    std::vector<char> stream_buffer{std::vector<char>(1'000'000)};

    /*!\name Stream / file access
     * \{
     */
    //!\brief The type of the internal stream pointers. Allows dynamically setting ownership management.
    using stream_ptr_t = std::unique_ptr<std::basic_ostream<stream_char_type>,
                                         std::function<void(std::basic_ostream<stream_char_type> *)>>;
    //!\brief Stream deleter that does nothing (no ownership assumed).
    static void stream_deleter_noop(std::basic_ostream<stream_char_type> *)
    {}
    //!\brief Stream deleter with default behaviour (ownership assumed).
    static void stream_deleter_default(std::basic_ostream<stream_char_type> * ptr)
    {
        delete ptr;
    }

    //!\brief The primary stream is the user provided stream or the file stream if constructed from filename.
    stream_ptr_t primary_stream{nullptr, stream_deleter_noop};
    //!\brief The secondary stream is a compression layer on the primary or just points to the primary (no compression).
    stream_ptr_t secondary_stream{nullptr, stream_deleter_noop};

    //!\brief Type of the format, a std::variant over the `valid_formats`.
    using format_type =
        typename detail::variant_from_tags<valid_formats, detail::sequence_file_output_format_exposer>::type;
    //!\brief The actual std::variant holding a pointer to the detected/selected format.
    format_type format;
    //!\}

    //!\brief Write record to format.
    template <typename seq_t, typename id_t, typename qual_t>
    void write_record(seq_t && seq, id_t && id, qual_t && qual)
    {
        assert(!format.valueless_by_exception());
        std::visit(
            [&](auto & f)
            {
                {
                    f.write_sequence_record(*secondary_stream, options, seq, id, qual);
                }
            },
            format);
    }

    //!\brief Befriend iterator so it can access the buffers.
    friend iterator;
};

/*!\name Type deduction guides
 * \relates seqan3::sequence_file_output
 * \{
 */

//!\brief Deduction guide for given stream and file format.
template <output_stream stream_t, sequence_file_output_format file_format>
sequence_file_output(stream_t &,
                     file_format const &)
    -> sequence_file_output<typename sequence_file_output<>::selected_field_ids, // default field ids
                            type_list<file_format>>;

//!\overload
template <output_stream stream_t, sequence_file_output_format file_format>
sequence_file_output(stream_t &&,
                     file_format const &)
    -> sequence_file_output<typename sequence_file_output<>::selected_field_ids, // default field ids.
                            type_list<file_format>>;

//!\brief Deduction guide for given stream, file format and field ids.
template <output_stream stream_t,
          sequence_file_output_format file_format,
          detail::fields_specialisation selected_field_ids>
sequence_file_output(stream_t &&,
                     file_format const &,
                     selected_field_ids const &) -> sequence_file_output<selected_field_ids, type_list<file_format>>;

//!\overload
template <output_stream stream_t,
          sequence_file_output_format file_format,
          detail::fields_specialisation selected_field_ids>
sequence_file_output(stream_t &,
                     file_format const &,
                     selected_field_ids const &) -> sequence_file_output<selected_field_ids, type_list<file_format>>;
//!\}
} // namespace seqan3

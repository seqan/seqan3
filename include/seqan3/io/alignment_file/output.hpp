// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::alignment_file_output and corresponding traits classes.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <cassert>
#include <fstream>
#include <string>
#include <string_view>
#include <variant>
#include <vector>

#include <seqan3/core/concept/tuple.hpp>
#include <seqan3/core/type_list/traits.hpp>
#include <seqan3/io/alignment_file/format_bam.hpp>
#include <seqan3/io/alignment_file/format_sam.hpp>
#include <seqan3/io/alignment_file/header.hpp>
#include <seqan3/io/alignment_file/misc.hpp>
#include <seqan3/io/alignment_file/output_format_concept.hpp>
#include <seqan3/io/alignment_file/output_options.hpp>
#include <seqan3/io/detail/out_file_iterator.hpp>
#include <seqan3/io/detail/misc_output.hpp>
#include <seqan3/io/detail/record.hpp>
#include <seqan3/io/exception.hpp>
#include <seqan3/std/filesystem>
#include <seqan3/io/record.hpp>
#include <seqan3/io/stream/concept.hpp>
#include <seqan3/std/ranges>

namespace seqan3
{

// ----------------------------------------------------------------------------
// alignment_file_output
// ----------------------------------------------------------------------------

/*!\brief A class for writing alignment files, e.g. SAM, BAL, BLAST, ...
 * \ingroup alignment_file
 * \tparam selected_field_ids   A seqan3::fields type with the list and order of
 *                              fields IDs; only relevant if these can't be deduced.
 * \tparam valid_formats        A seqan3::type_list of the selectable formats (each
 *                              must model seqan3::alignment_file_output_format).
 *
 * \details
 *
 * \copydetails alignment_file
 *
 * The member functions take any and either of these fields.
 *
 * ### Construction and specialisation
 *
 * This class comes with two constructors, one for construction from a file name
 * and one for construction from an existing stream and a known format. The first
 * one automatically picks the format based on the extension of the file name.
 * The second can be used if you have a non-file stream, like std::cout or
 * std::ostringstream, that you want to read from and/or if you cannot use
 * file-extension based detection, but know that your output file has a certain
 * format.
 *
 * In most cases the template parameters are deduced completely automatically:
 *
 * \include test/snippet/io/alignment_file/alignment_file_output_filename_construction.cpp
 *
 * Writing to std::cout:
 *
 * \include test/snippet/io/alignment_file/alignment_file_output_cout_write.cpp
 *
 * Note that this is not the same as writing `alignment_file_output<>`
 * (with angle brackets). In the latter case they are explicitly set to their
 * default values, in the former case
 * [automatic deduction](https://en.cppreference.com/w/cpp/language/class_template_argument_deduction)
 * happens which chooses different parameters depending on the constructor arguments.
 * For opening from file, `alignment_file_output<>` would have also worked, but for
 * opening from stream it would not have.
 *
 * ### Writing record-wise
 *
 * \include test/snippet/io/alignment_file/record_based_writing.cpp
 *
 * The easiest way to write to an alignment file is to use the push_back() member functions. These
 * work similarly to how they work on an std::vector.
 * You may also use a tuple like interface or the emplace_back()
 * function but this is not recommended since one would have to keep track of the
 * correct order of many fields (14 in total). For the record based interface
 * using push_back please also see the seqan3::record documentation on how to specify
 * a record with the correct field and type lists.
 *
 * You may also use the output file's iterator for writing, however, this rarely provides an advantage.
 *
 * ### Writing record-wise (custom fields)
 *
 * If you want to omit non-required parameter or
 * change the order of the parameters, you can pass a non-empty fields trait object to the
 * seqan3::alignment_file_output constructor to select the fields that are used for interpreting the arguments.
 *
 * The following snippets demonstrates the usage of such a field_traits object.
 *
 * \include test/snippet/io/alignment_file/record_based_writing2.cpp
 *
 * A different way of passing custom fields to the file is to pass a seqan3::record – instead of a tuple – to
 * push_back(). The seqan3::record clearly indicates which of its elements has which seqan3::field so **the file will
 * use that information instead of the template argument**. This is especially handy when reading from one file and
 * writing to another, because you don't have to configure the output file to match the input file, it will just work:
 *
 * \include test/snippet/io/alignment_file/alignment_file_output_custom_fields.cpp
 *
 * This will copy the FLAG and REF_OFFSET value into the new output file. Note that the other SAM columns in the
 * output file will be defaulted, so unless you specify to read all SAM columns (see seqan3::format_sam)
 * the output file will not be equal to the input file.
 *
 * ### Writing record-wise in batches
 *
 * You can write multiple records at once, by assigning to the file:
 *
 * \include test/snippet/io/alignment_file/alignment_file_output_write_range.cpp
 *
 * ### File I/O pipelines
 *
 * Record-wise writing in batches also works for writing from input files directly to output files, because input
 * files are also input ranges in SeqAn:
 *
 * \include test/snippet/io/alignment_file/alignment_file_output_input_range.cpp
 *
 * This can be combined with file-based views to create I/O pipelines:
 *
 * \include test/snippet/io/alignment_file/alignment_file_output_io_pipeline.cpp
 *
 * ### Formats
 *
 * We currently support writing the following formats:
 *   * seqan3::format_sam
 *   * seqan3::format_bam
 */
template <detail::fields_specialisation selected_field_ids_ =
              fields<field::seq,
                     field::id,
                     field::offset,
                     field::ref_seq,
                     field::ref_id,
                     field::ref_offset,
                     field::alignment,
                     field::cigar,
                     field::mapq,
                     field::qual,
                     field::flag,
                     field::mate,
                     field::tags,
                     field::evalue,
                     field::bit_score,
                     field::header_ptr>,
          detail::type_list_of_alignment_file_output_formats valid_formats_ = type_list<format_sam, format_bam>,
          typename ref_ids_type = ref_info_not_given>
class alignment_file_output
{
public:
    /*!\name Template arguments
     * \brief Exposed as member types for public access.
     * \{
     */
    //!\brief A seqan3::fields list with the fields selected for the record.
    using selected_field_ids    = selected_field_ids_;
    //!\brief A seqan3::type_list with the possible formats.
    using valid_formats         = valid_formats_;
    //!\brief Character type of the stream(s).
    using stream_char_type      = char;
    //!\}

    //!\brief The subset of seqan3::field IDs that are valid for this file.
    using field_ids             = fields<field::header_ptr,
                                         field::seq,
                                         field::id,
                                         field::offset,
                                         field::ref_seq,
                                         field::ref_id,
                                         field::ref_offset,
                                         field::alignment,
                                         field::cigar,
                                         field::mapq,
                                         field::flag,
                                         field::qual,
                                         field::mate,
                                         field::tags,
                                         field::evalue,
                                         field::bit_score>;

    static_assert([] () constexpr
                  {
                      for (field f : selected_field_ids::as_array)
                          if (!field_ids::contains(f))
                              return false;
                      return true;
                  }(),
                  "You selected a field that is not valid for alignment files, "
                  "please refer to the documentation of "
                  "seqan3::alignment_file_output::field_ids for the accepted values.");

    /*!\name Range associated types
     * \brief Most of the range associated types are `void` for output ranges.
     * \{
     */

    //!\brief The value type (void).
    using value_type        = void;
    //!\brief The reference type (void).
    using reference         = void;
    //!\brief The const reference type (void).
    using const_reference   = void;
    //!\brief The size type (void).
    using size_type         = void;
    //!\brief A signed integer type, usually std::ptrdiff_t.
    using difference_type   = std::ptrdiff_t;
    //!\brief The iterator type of this view (an output iterator).
    using iterator          = detail::out_file_iterator<alignment_file_output>;
    //!\brief The const iterator type is void, because files are not const-iterable.
    using const_iterator    = void;
    //!\brief The type returned by end().
    using sentinel          = std::ranges::default_sentinel_t;
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    //!\brief Default constructor is explicitly deleted, you need to give a stream or file name.
    alignment_file_output() = delete;
    //!\brief Copy construction is explicitly deleted, because you can't have multiple access to the same file.
    alignment_file_output(alignment_file_output const &) = delete;
    //!\brief Copy assignment is explicitly deleted, because you can't have multiple access to the same file.
    alignment_file_output & operator=(alignment_file_output const &) = delete;
    //!\brief Move construction is defaulted.
    alignment_file_output(alignment_file_output &&) = default;
    //!\brief Move assignment is defaulted.
    alignment_file_output & operator=(alignment_file_output &&) = default;
    //!\brief Destructor is defaulted.
    ~alignment_file_output() = default;

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
     *
     * ### Example:
     *
     * In most cases the template parameters are deduced completely automatically:
     *
     * \include test/snippet/io/alignment_file/alignment_file_output_filename_construction.cpp
     *
     * Writing with custom selected fields:
     *
     * \include test/snippet/io/alignment_file/alignment_file_output_format_construction.cpp
     */
    alignment_file_output(std::filesystem::path filename,
                          selected_field_ids const & SEQAN3_DOXYGEN_ONLY(fields_tag) = selected_field_ids{}) :
        primary_stream{new std::ofstream{}, stream_deleter_default}
    {
        primary_stream->rdbuf()->pubsetbuf(stream_buffer.data(), stream_buffer.size());
        static_cast<std::basic_ofstream<char> *>(primary_stream.get())->open(filename,
                                                                             std::ios_base::out | std::ios::binary);

        // open stream
        if (!primary_stream->good())
            throw file_open_error{"Could not open file " + filename.string() + " for writing."};

        // possibly add intermediate compression stream
        secondary_stream = detail::make_secondary_ostream(*primary_stream, filename);

        // initialise format handler or throw if format is not found
        detail::set_format(format, filename);
    }

    /*!\brief Construct from an existing stream and with specified format.
     * \tparam stream_type   The type of stream to write to; must model seqan3::output_stream.
     * \tparam file_format   The format of the file in the stream, must satisfy seqan3::alignment_file_output_format.
     * \param[out] stream    The stream to write to, must be derived of std::basic_ostream<stream_char_t>.
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
    template <output_stream stream_type, alignment_file_output_format file_format>
    //!\cond
        requires std::same_as<typename std::remove_reference_t<stream_type>::char_type, stream_char_type>
    //!\endcond
    alignment_file_output(stream_type              & stream,
                          file_format        const & SEQAN3_DOXYGEN_ONLY(format_tag),
                          selected_field_ids const & SEQAN3_DOXYGEN_ONLY(fields_tag) = selected_field_ids{}) :
        primary_stream{&stream, stream_deleter_noop},
        secondary_stream{&stream, stream_deleter_noop},
        format{detail::alignment_file_output_format_exposer<file_format>{}}
    {
        static_assert(list_traits::contains<file_format, valid_formats>,
                      "You selected a format that is not in the valid_formats of this file.");
    }

    //!\overload
    template <output_stream stream_type, alignment_file_output_format file_format>
    //!\cond
        requires std::same_as<typename std::remove_reference_t<stream_type>::char_type, stream_char_type>
    //!\endcond
    alignment_file_output(stream_type             && stream,
                          file_format        const & SEQAN3_DOXYGEN_ONLY(format_tag),
                          selected_field_ids const & SEQAN3_DOXYGEN_ONLY(fields_tag) = selected_field_ids{}) :
        primary_stream{new stream_type{std::move(stream)}, stream_deleter_default},
        secondary_stream{&*primary_stream, stream_deleter_noop},
        format{detail::alignment_file_output_format_exposer<file_format>{}}
    {
        static_assert(list_traits::contains<file_format, valid_formats>,
                      "You selected a format that is not in the valid_formats of this file.");
    }

    /*!\brief Construct from filename.
     * \tparam ref_ids_type_    The type of range over reference ids; must model std::forward_range.
     * \tparam ref_lengths_type The type of range over reference lengths; must model std::forward_range.
     *
     * \param[in] filename      Path to the file you wish to open.
     * \param[in] ref_ids       A range over reference ids.
     * \param[in] ref_lengths   A range over lengths of reference sequences (same order as ref_ids).
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
     *
     * ### Example:
     *
     * In most cases the template parameters are deduced completely automatically:
     *
     * \include test/snippet/io/alignment_file/alignment_file_output_filename_construction_with_ref_info.cpp
     *
     * Writing with custom selected fields:
     *
     * \include test/snippet/io/alignment_file/alignment_file_output_format_construction.cpp
     */
    template <typename ref_ids_type_, std::ranges::forward_range ref_lengths_type>
    //!\cond
        requires std::same_as<std::remove_reference_t<ref_ids_type_>, ref_ids_type>
    //!\endcond
    alignment_file_output(std::filesystem::path const & filename,
                          ref_ids_type_              && ref_ids,
                          ref_lengths_type           && ref_lengths,
                          selected_field_ids    const & SEQAN3_DOXYGEN_ONLY(fields_tag) = selected_field_ids{}) :
        alignment_file_output{filename, selected_field_ids{}}

    {
        initialise_header_information(ref_ids, ref_lengths);
    }

    /*!\brief Construct from an existing stream and with specified format.
     * \tparam stream_type      The type of stream to write to; must model seqan3::output_stream.
     * \tparam file_format      The format of the file in the stream, must model seqan3::alignment_file_output_format.
     * \tparam ref_ids_type_    The type of range over reference ids; must model std::forward_range.
     * \tparam ref_lengths_type The type of range over reference lengths; must model std::forward_range.
     *
     * \param[in] stream        The stream to operate on (this must be std::move'd in!).
     * \param[in] ref_ids       A range over reference ids.
     * \param[in] ref_lengths   A range over lengths of reference sequences (same order as ref_ids).
     * \param[in] format_tag    The file format tag.
     * \param[in] fields_tag    A seqan3::fields tag. [optional]
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
    template <output_stream stream_type,
              alignment_file_output_format file_format,
              typename ref_ids_type_, // generic type to capture lvalue references
              std::ranges::forward_range ref_lengths_type>
    //!\cond
        requires std::same_as<std::remove_reference_t<ref_ids_type_>, ref_ids_type>
    //!\endcond
    alignment_file_output(stream_type             && stream,
                          ref_ids_type_           && ref_ids,
                          ref_lengths_type        && ref_lengths,
                          file_format        const & SEQAN3_DOXYGEN_ONLY(format_tag),
                          selected_field_ids const & SEQAN3_DOXYGEN_ONLY(fields_tag) = selected_field_ids{}) :
        alignment_file_output{std::forward<stream_type>(stream), file_format{}, selected_field_ids{}}
    {
        initialise_header_information(ref_ids, ref_lengths);
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
     * \include test/snippet/io/alignment_file/begin_iterator.cpp
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
     * Constant.
     *
     * ### Exceptions
     *
     * Basic exception safety.
     *
     * ### Example
     *
     * \include test/snippet/io/alignment_file/push_back_record.cpp
     */
    template <typename record_t>
    void push_back(record_t && r)
    //!\cond
        requires tuple_like<record_t> &&
                 requires { requires detail::is_type_specialisation_of_v<remove_cvref_t<record_t>, record>; }
    //!\endcond
    {
        using default_align_t = std::pair<std::span<gapped<char>>, std::span<gapped<char>>>;
        using default_mate_t  = std::tuple<std::string_view, std::optional<int32_t>, int32_t>;

        write_record(detail::get_or<field::header_ptr>(r, nullptr),
                     detail::get_or<field::seq>(r, std::string_view{}),
                     detail::get_or<field::qual>(r, std::string_view{}),
                     detail::get_or<field::id>(r, std::string_view{}),
                     detail::get_or<field::offset>(r, 0u),
                     detail::get_or<field::ref_seq>(r, std::string_view{}),
                     detail::get_or<field::ref_id>(r, std::ignore),
                     detail::get_or<field::ref_offset>(r, std::optional<int32_t>{}),
                     detail::get_or<field::alignment>(r, default_align_t{}),
                     detail::get_or<field::cigar>(r, std::vector<cigar>{}),
                     detail::get_or<field::flag>(r, sam_flag::none),
                     detail::get_or<field::mapq>(r, 0u),
                     detail::get_or<field::mate>(r, default_mate_t{}),
                     detail::get_or<field::tags>(r, sam_tag_dictionary{}),
                     detail::get_or<field::evalue>(r, 0u),
                     detail::get_or<field::bit_score>(r, 0u));
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
     * Constant.
     *
     * ### Exceptions
     *
     * Basic exception safety.
     *
     * ### Example
     *
     * \include test/snippet/io/alignment_file/push_back_tuple.cpp
     */
    template <typename tuple_t>
    void push_back(tuple_t && t)
    //!\cond
        requires tuple_like<tuple_t>
    //!\endcond
    {
        using default_align_t = std::pair<std::span<gapped<char>>, std::span<gapped<char>>>;
        using default_mate_t  = std::tuple<std::string_view, std::optional<int32_t>, int32_t>;

        // index_of might return npos, but this will be handled well by get_or_ignore (and just return ignore)
        write_record(detail::get_or<selected_field_ids::index_of(field::header_ptr)>(t, nullptr),
                     detail::get_or<selected_field_ids::index_of(field::seq)>(t, std::string_view{}),
                     detail::get_or<selected_field_ids::index_of(field::qual)>(t, std::string_view{}),
                     detail::get_or<selected_field_ids::index_of(field::id)>(t, std::string_view{}),
                     detail::get_or<selected_field_ids::index_of(field::offset)>(t, 0u),
                     detail::get_or<selected_field_ids::index_of(field::ref_seq)>(t, std::string_view{}),
                     detail::get_or<selected_field_ids::index_of(field::ref_id)>(t, std::ignore),
                     detail::get_or<selected_field_ids::index_of(field::ref_offset)>(t, std::optional<int32_t>{}),
                     detail::get_or<selected_field_ids::index_of(field::alignment)>(t, default_align_t{}),
                     detail::get_or<selected_field_ids::index_of(field::cigar)>(t, std::vector<cigar>{}),
                     detail::get_or<selected_field_ids::index_of(field::flag)>(t, sam_flag::none),
                     detail::get_or<selected_field_ids::index_of(field::mapq)>(t, 0u),
                     detail::get_or<selected_field_ids::index_of(field::mate)>(t, default_mate_t{}),
                     detail::get_or<selected_field_ids::index_of(field::tags)>(t, sam_tag_dictionary{}),
                     detail::get_or<selected_field_ids::index_of(field::evalue)>(t, 0u),
                     detail::get_or<selected_field_ids::index_of(field::bit_score)>(t, 0u));
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
     * Constant.
     *
     * ### Exceptions
     *
     * Basic exception safety.
     *
     * ### Example
     *
     * \include test/snippet/io/alignment_file/emplace_back.cpp
     */
    template <typename arg_t, typename ...arg_types>
    void emplace_back(arg_t && arg, arg_types && ... args)
    {
        push_back(std::tie(arg, args...));
    }

    /*!\brief            Write a range of records (or tuples) to the file.
     * \tparam rng_t     Type of the range, must satisfy seqan3::output_range and have a reference type that
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
     * \include test/snippet/io/alignment_file/alignment_file_output_write_range.cpp
     */
    template <typename rng_t>
    alignment_file_output & operator=(rng_t && range)
    //!\cond
        requires std::ranges::input_range<rng_t> && tuple_like<std::ranges::range_reference_t<rng_t>>
    //!\endcond
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
     * This operator enables alignment_file_output to be at the end of a piping operation. It just calls
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
     * \include test/snippet/io/alignment_file/alignment_file_output_write_range.cpp
     *
     * This is especially useful in combination with file-based filters:
     *
     * \include test/snippet/io/alignment_file/alignment_file_output_io_pipeline.cpp
     *
     */
    template <typename rng_t>
    friend alignment_file_output & operator|(rng_t && range, alignment_file_output & f)
    //!\cond
        requires std::ranges::input_range<rng_t> && tuple_like<std::ranges::range_reference_t<rng_t>>
    //!\endcond
    {
        f = range;
        return f;
    }

    //!\overload
    template <typename rng_t>
    friend alignment_file_output operator|(rng_t && range, alignment_file_output && f)
    //!\cond
        requires std::ranges::input_range<rng_t> && tuple_like<std::ranges::range_reference_t<rng_t>>
    //!\endcond
    {
        f = range;
        return std::move(f);
    }
    //!\}

    //!\brief The options are public and its members can be set directly.
    alignment_file_output_options options;

    /*!\cond DEV
     * \brief Expose a reference to the secondary stream. [public, but not documented as part of the API]
     */
    std::basic_ostream<stream_char_type> & get_stream()
    {
        return *secondary_stream;
    }
    //!\endcond

    /*!\brief Access the file's header.
     *
     * \details
     *
     * ### Example
     *
     * \include test/snippet/io/alignment_file/alignment_file_output_set_header.cpp
     *
     * \sa seqan3::alignment_file_header
     */
    auto & header()
    {
        if constexpr (std::same_as<ref_ids_type, ref_info_not_given>)
            throw std::logic_error{"Please construct your file with reference id and length information in order "
                                   "to properly initialise the header before accessing it."};

        return *header_ptr;
    }

protected:
    //!\privatesection
    //!\brief A larger (compared to stl default) stream buffer to use when reading from a file.
    std::vector<char> stream_buffer{std::vector<char>(1'000'000)};

    /*!\name Stream / file access
     * \{
     */
    //!\brief The type of the internal stream pointers. Allows dynamically setting ownership management.
    using stream_ptr_t = std::unique_ptr<std::basic_ostream<stream_char_type>,
                                         std::function<void(std::basic_ostream<stream_char_type>*)>>;
    //!\brief Stream deleter that does nothing (no ownership assumed).
    static void stream_deleter_noop(std::basic_ostream<stream_char_type> *) {}
    //!\brief Stream deleter with default behaviour (ownership assumed).
    static void stream_deleter_default(std::basic_ostream<stream_char_type> * ptr) { delete ptr; }

    //!\brief The primary stream is the user provided stream or the file stream if constructed from filename.
    stream_ptr_t primary_stream{nullptr, stream_deleter_noop};
    //!\brief The secondary stream is a compression layer on the primary or just points to the primary (no compression).
    stream_ptr_t secondary_stream{nullptr, stream_deleter_noop};

    //!\brief Type of the format, an std::variant over the `valid_formats`.
    using format_type = typename detail::variant_from_tags<valid_formats,
                                                           detail::alignment_file_output_format_exposer>::type;

    //!\brief The actual std::variant holding a pointer to the detected/selected format.
    format_type format;
    //!\}

    //!\brief The header type, which specilised with ref_ids_type if reference information are given.
    using header_type = alignment_file_header<std::conditional_t<std::same_as<ref_ids_type, ref_info_not_given>,
                                              std::vector<std::string>,
                                              ref_ids_type>>;

    //!\brief The file header object (will be set on construction).
    std::unique_ptr<header_type> header_ptr;

    //!\brief Fill the header reference dictionary, with the given info.
    template <typename ref_ids_type_, typename ref_lengths_type>
    void initialise_header_information(ref_ids_type_ && ref_ids, ref_lengths_type && ref_lengths)
    {
        assert(std::ranges::size(ref_ids) == std::ranges::size(ref_lengths));

        header_ptr = std::make_unique<alignment_file_header<ref_ids_type>>(std::forward<ref_ids_type_>(ref_ids));

        for (int32_t idx = 0; idx < std::ranges::distance(ref_ids); ++idx)
        {
            header_ptr->ref_id_info.emplace_back(ref_lengths[idx], "");

            if constexpr (std::ranges::contiguous_range<std::ranges::range_reference_t<ref_ids_type_>> &&
                          std::ranges::sized_range<std::ranges::range_reference_t<ref_ids_type_>> &&
                          forwarding_range<std::ranges::range_reference_t<ref_ids_type_>>)
            {
                auto && id = header_ptr->ref_ids()[idx];
                header_ptr->ref_dict[std::span{std::ranges::data(id), std::ranges::size(id)}] = idx;
            }
            else
            {
                header_ptr->ref_dict[header_ptr->ref_ids()[idx]] = idx;
            }
        }
    }

    //!\brief Write record to format.
    template <typename record_header_ptr_t, typename ...pack_type>
    void write_record(record_header_ptr_t && record_header_ptr, pack_type && ...remainder)
    {
        static_assert((sizeof...(pack_type) == 15), "Wrong parameter list passed to write_record.");

        assert(!format.valueless_by_exception());

        std::visit([&] (auto & f)
        {
            // use header from record if explicitly given, e.g. file_output = file_input
            if constexpr (!std::same_as<record_header_ptr_t, std::nullptr_t>)
            {
                f.write_alignment_record(*secondary_stream,
                                         options,
                                         *record_header_ptr,
                                         std::forward<pack_type>(remainder)...);
            }
            else if constexpr (std::same_as<ref_ids_type, ref_info_not_given>)
            {
                f.write_alignment_record(*secondary_stream,
                                         options,
                                         std::ignore,
                                         std::forward<pack_type>(remainder)...);
            }
            else
            {
                f.write_alignment_record(*secondary_stream,
                                         options,
                                         *header_ptr,
                                         std::forward<pack_type>(remainder)...);
            }
        }, format);
    }

    //!\brief Befriend iterator so it can access the buffers.
    friend iterator;
};

/*!\name Type deduction guides
 * \relates seqan3::alignment_file_output
 * \{
 */

/*!\brief Deduces selected_field_ids from input and sets alignment_file_output::ref_ids_type to
 *        seqan3::detail::ref_info_not_given. valid_formats is defaulted.
 */
template <detail::fields_specialisation selected_field_ids>
alignment_file_output(std::filesystem::path, selected_field_ids const &)
    -> alignment_file_output<selected_field_ids,
                             typename alignment_file_output<>::valid_formats,
                             ref_info_not_given>;

/*!\brief Deduces selected_field_ids, and the valid format from input and
 *        sets alignment_file_output::ref_ids_type to seqan3::detail::ref_info_not_given.
 */
template <output_stream stream_type,
          alignment_file_output_format file_format,
          detail::fields_specialisation selected_field_ids>
alignment_file_output(stream_type &&, file_format const &, selected_field_ids const &)
    -> alignment_file_output<selected_field_ids,
                             type_list<file_format>,
                             ref_info_not_given>;

/*!\brief Deduces selected_field_ids, and the valid format from input and
 *        sets alignment_file_output::ref_ids_type to seqan3::detail::ref_info_not_given.
 */
template <output_stream stream_type,
          alignment_file_output_format file_format,
          detail::fields_specialisation selected_field_ids>
alignment_file_output(stream_type &, file_format const &, selected_field_ids const &)
    -> alignment_file_output<selected_field_ids,
                             type_list<file_format>,
                             ref_info_not_given>;

/*!\brief Deduces the valid format from input and sets alignment_file_output::ref_ids_type to
 *        seqan3::detail::ref_info_not_given. selected_field_ids is defaulted.
 */
template <output_stream stream_type,
          alignment_file_output_format file_format>
alignment_file_output(stream_type &&, file_format const &)
    -> alignment_file_output<typename alignment_file_output<>::selected_field_ids,
                             type_list<file_format>,
                             ref_info_not_given>;

/*!\brief Deduces the valid format from input and sets alignment_file_output::ref_ids_type to
 *        seqan3::detail::ref_info_not_given. selected_field_ids is defaulted.
 */
template <output_stream stream_type,
          alignment_file_output_format file_format>
alignment_file_output(stream_type &, file_format const &)
    -> alignment_file_output<typename alignment_file_output<>::selected_field_ids,
                             type_list<file_format>,
                             ref_info_not_given>;

//!\brief Deduces selected_field_ids and ref_ids_type from input. valid_formats is defaulted.
template <detail::fields_specialisation selected_field_ids,
          std::ranges::forward_range ref_ids_type,
          std::ranges::forward_range ref_lengths_type>
alignment_file_output(std::filesystem::path const &,
                      ref_ids_type &&,
                      ref_lengths_type &&,
                      selected_field_ids const &)
    -> alignment_file_output<selected_field_ids,
                             typename alignment_file_output<>::valid_formats,
                             std::remove_reference_t<ref_ids_type>>;

//!\brief Deduces ref_ids_type from input. Valid formats, and selected_field_ids are defaulted.
template <std::ranges::forward_range ref_ids_type,
          std::ranges::forward_range ref_lengths_type>
alignment_file_output(std::filesystem::path const &,
                      ref_ids_type &&,
                      ref_lengths_type &&)
    -> alignment_file_output<typename alignment_file_output<>::selected_field_ids,
                             typename alignment_file_output<>::valid_formats,
                             std::remove_reference_t<ref_ids_type>>;

//!\brief Deduces selected_field_ids, the valid format, and the ref_ids_type from input.
template <output_stream stream_type,
          std::ranges::forward_range ref_ids_type,
          std::ranges::forward_range ref_lengths_type,
          alignment_file_output_format file_format,
          detail::fields_specialisation selected_field_ids>
alignment_file_output(stream_type &&,
                      ref_ids_type &&,
                      ref_lengths_type &&,
                      file_format const &,
                      selected_field_ids const &)
    -> alignment_file_output<selected_field_ids,
                             type_list<file_format>,
                             std::remove_reference_t<ref_ids_type>>;

//!\brief Deduces selected_field_ids, the valid format, and the ref_ids_type from input.
template <output_stream stream_type,
          std::ranges::forward_range ref_ids_type,
          std::ranges::forward_range ref_lengths_type,
          alignment_file_output_format file_format,
          detail::fields_specialisation selected_field_ids>
alignment_file_output(stream_type &,
                      ref_ids_type &&,
                      ref_lengths_type &&,
                      file_format const &,
                      selected_field_ids const &)
    -> alignment_file_output<selected_field_ids,
                             type_list<file_format>,
                             std::remove_reference_t<ref_ids_type>>;

//!\brief Deduces the valid format, and the ref_ids_type from input. selected_field_ids is defaulted.
template <output_stream stream_type,
          std::ranges::forward_range ref_ids_type,
          std::ranges::forward_range ref_lengths_type,
          alignment_file_output_format file_format>
alignment_file_output(stream_type &&,
                      ref_ids_type &&,
                      ref_lengths_type &&,
                      file_format const &)
    -> alignment_file_output<typename alignment_file_output<>::selected_field_ids,
                             type_list<file_format>,
                             std::remove_reference_t<ref_ids_type>>;

//!\brief Deduces the valid format, and the ref_ids_type from input. selected_field_ids is defaulted.
template <output_stream stream_type,
          std::ranges::forward_range ref_ids_type,
          std::ranges::forward_range ref_lengths_type,
          alignment_file_output_format file_format>
alignment_file_output(stream_type &,
                      ref_ids_type &&,
                      ref_lengths_type &&,
                      file_format const &)
    -> alignment_file_output<typename alignment_file_output<>::selected_field_ids,
                             type_list<file_format>,
                             std::remove_reference_t<ref_ids_type>>;
//!\}

} // namespace seqan3

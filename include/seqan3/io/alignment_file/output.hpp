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
 * \brief Provides seqan3::alignment_file_output and corresponding traits classes.
 * \author Svenja Mehringer <avenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <cassert>
#include <fstream>
#include <string>
#include <string_view>
#include <variant>
#include <vector>

#include <seqan3/core/concept/tuple.hpp>
#include <seqan3/core/metafunction/basic.hpp>
#include <seqan3/core/metafunction/template_inspection.hpp>
#include <seqan3/io/alignment_file/format_sam.hpp>
#include <seqan3/io/alignment_file/header.hpp>
#include <seqan3/io/alignment_file/OutputFormat.hpp>
#include <seqan3/io/alignment_file/output_options.hpp>
#include <seqan3/io/detail/out_file_iterator.hpp>
#include <seqan3/io/detail/record.hpp>
#include <seqan3/io/exception.hpp>
#include <seqan3/io/filesystem.hpp>
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
 *                              must model seqan3::AlignmentFileOutputFormat).
 * \tparam stream_type          The type of the stream, must model seqan3::Ostream.
 *
 * \details
 *
 * ### Introduction
 *
 * Alignment files provide a way to store pairwise alignment information.
 * Well-known formats include SAM and BAM or BLAST.
 *
 * The Alignment file abstraction supports writing following fields:
 *
 *   1. field::SEQ
 *   2. field::ID
 *   3. field::OFFSET
 *   4. field::REF_SEQ
 *   5. field::REF_ID
 *   6. field::REF_OFFSET
 *   7. field::ALIGNMENT
 *   8. field::MAPQ
 *   9. field::FLAG
 *   10. field::QUAL
 *   11. field::MATE
 *   12. field::TAGS
 *   13. field::EVALUE
 *   14. field::BIT_SCORE
 *
 * There is an additional field called seqan3::field::HEADER_PTR. It is used to transfer
 * header information from seqan3::alignment_file_input to seqan3::alignment_file_output,
 * but you needn't deal with this field manually.
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
 * \snippet test/snippet/io/alignment_file/output.cpp filename_construction
 *
 * Writing to std::cout:
 *
 * \snippet test/snippet/io/alignment_file/output.cpp format_construction
 *
 * Note that this is not the same as writing `alignment_file_output<>`
 * (with angle brackets). In the latter case they are explicitly set to their
 * default values, in the former case
 * [automatic deduction](http://en.cppreference.com/w/cpp/language/class_template_argument_deduction)
 * happens which chooses different parameters depending on the constructor arguments.
 * For opening from file, `alignment_file_output<>` would have also worked, but for
 * opening from stream it would not have.
 *
 * ### Writing record-wise
 *
 * \snippet test/snippet/io/alignment_file/record_based_writing.cpp all
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
 * \snippet test/snippet/io/alignment_file/record_based_writing2.cpp all
 *
 * A different way of passing custom fields to the file is to pass a seqan3::record – instead of a tuple – to
 * push_back(). The seqan3::record clearly indicates which of its elements has which seqan3::field so **the file will
 * use that information instead of the template argument**. This is especially handy when reading from one file and
 * writing to another, because you don't have to configure the output file to match the input file, it will just work:
 *
 * \todo TODO introduce snippet once seqan3:;alignment_file_in is implemented.
 *
 * ```cpp
 * alignment_file_in  fin{"input.sam", fields<field::REF_OFFSET, field::FLAG>{}};
 * alignment_file_output fout{"output.sam"}; // doesn't have to match the configuration
 *
 * for (auto & r : fin)
 * {
 *         fout.push_back(r); // copy all the records.
 * }
 * ```
 * This will copy the FLAG and REF_OFFSET value into the new output file. Note that the other SAM columns in the
 * output file will be defaulted, so unless you specify to read all SAM columns (see seqan3::alignment_file_format_sam)
 * the output file will not be equal to the input file.
 *
 * ### Writing record-wise in batches
 *
 * You can write multiple records at once, by assigning to the file:
 *
 * \snippet test/snippet/io/alignment_file/output.cpp write_range
 *
 * ### File I/O pipelines
 *
 * Record-wise writing in batches also works for writing from input files directly to output files, because input
 * files are also input ranges in SeqAn3:
 *
 * \todo TODO introduce snippet once seqan3:;alignment_file_in is implemented.
 *
 * ```cpp
 * // file format conversion in one line:
 * alignment_file_output fout{"output.sam"} = sequence_file_in{"input.sam"};
 *
 * // or in pipe notation:
 * sequence_file_in{"input.sam"} | alignment_file_output{"output.sam"};
 * ```
 *
 * This can be combined with file-based views to create I/O pipelines:
 *
 * ```cpp
 * sequence_file_in{"input.sam"} | ranges::view::take(5) // take only the first 5 records
 *                               | alignment_file_output{"output.sam"};
 * ```
 *
 * ### Formats
 *
 * TODO give overview of formats, once they are all implemented
 *
 * \sa seqan3::alignment_file_format_sam
 */
template <detail::Fields  selected_field_ids_ =
              fields<field::SEQ,
                     field::ID,
                     field::OFFSET,
                     field::REF_SEQ,
                     field::REF_ID,
                     field::REF_OFFSET,
                     field::ALIGNMENT,
                     field::MAPQ,
                     field::QUAL,
                     field::FLAG,
                     field::MATE,
                     field::TAGS,
                     field::EVALUE,
                     field::BIT_SCORE,
                     field::HEADER_PTR>,
          detail::TypeListOfAlignmentFileOutputFormats valid_formats_ =
              type_list<alignment_file_format_sam/*,
                        alignment_file_format_bam,
                        alignment_file_format_blast_tabular*/>,
          Ostream<char> stream_type_ = std::ofstream>
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
    //!\brief The type of the underlying stream.
    using stream_type           = stream_type_;
    //!\}

    //!\brief The subset of seqan3::field IDs that are valid for this file.
    using field_ids             = fields<field::HEADER_PTR,
                                         field::SEQ,
                                         field::ID,
                                         field::OFFSET,
                                         field::REF_SEQ,
                                         field::REF_ID,
                                         field::REF_OFFSET,
                                         field::ALIGNMENT,
                                         field::MAPQ,
                                         field::FLAG,
                                         field::QUAL,
                                         field::MATE,
                                         field::TAGS,
                                         field::EVALUE,
                                         field::BIT_SCORE>;

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
    using value_type        = void;
    using reference         = void;
    using const_reference   = void;
    using size_type         = void;
    //!\brief A signed integer type, usually std::ptrdiff_t.
    using difference_type   = std::ptrdiff_t;
    //!\brief The iterator type of this view (an output iterator).
    using iterator          = detail::out_file_iterator<alignment_file_output>;
    //!\brief The const iterator type is void, because files are not const-iterable.
    using const_iterator    = void;
    //!\brief The type returned by end().
    using sentinel          = std::ranges::default_sentinel;
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
     * \param[in] _file_name    Path to the file you wish to open.
     * \param[in] fields_tag    A seqan3::fields tag. [optional]
     *
     * \details
     *
     * In addition to the file name, you may specify a custom seqan3::fields type which may be easier than
     * defining all the template parameters.
     *
     * ### Example:
     *
     * In most cases the template parameters are deduced completely automatically:
     *
     * \snippet test/snippet/io/alignment_file/output.cpp filename_construction
     *
     * Writing with custom selected fields:
     *
     * \snippet test/snippet/io/alignment_file/output.cpp format_construction
     */
    alignment_file_output(filesystem::path const & _file_name,
                          selected_field_ids const & SEQAN3_DOXYGEN_ONLY(fields_tag) = selected_field_ids{})
    {
        // open stream
        stream.open(_file_name, std::ios_base::out | std::ios::binary);
        if (!stream.is_open())
            throw file_open_error{"Could not open file for writing."};

        // initialise format handler or throw if format is not found
        detail::set_format(format, _file_name);
    }

    /*!\brief Construct from an existing stream and with specified format.
     * \tparam file_format   The format of the file in the stream, must satisfy seqan3::AlignmentFileOutputFormat.
     * \param[in] _stream    The stream to operate on (this must be std::move'd in!).
     * \param[in] format_tag The file format tag.
     * \param[in] fields_tag A seqan3::fields tag. [optional]
     *
     * \details
     */
    template <AlignmentFileOutputFormat file_format>
    alignment_file_output(stream_type             && _stream,
                          file_format        const & SEQAN3_DOXYGEN_ONLY(format_tag),
                          selected_field_ids const & SEQAN3_DOXYGEN_ONLY(fields_tag) = selected_field_ids{}) :
        stream{std::move(_stream)}, format{file_format{}}
    {
        static_assert(meta::in<valid_formats, file_format>::value,
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
     * \snippet test/snippet/io/alignment_file/begin_iterator.cpp all
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
     * Constant. TODO linear in the size of the written alignments?
     *
     * ### Exceptions
     *
     * Basic exception safety.
     *
     * ### Example
     *
     * \snippet test/snippet/io/alignment_file/push_back_record.cpp all
     */
    template <typename record_t>
    void push_back(record_t && r)
    //!\cond
        requires TupleLike<record_t> &&
                 requires { requires detail::is_type_specialisation_of_v<remove_cvref_t<record_t>, record>; }
    //!\endcond
    {
        using default_align_t = std::pair<std::basic_string_view<gapped<char>>, std::basic_string_view<gapped<char>>>;
        using default_mate_t  = std::tuple<std::string_view, uint32_t, uint32_t>;

        write_record(detail::get_or<field::HEADER_PTR>(r, nullptr),
                     detail::get_or<field::SEQ>(r, std::string_view{}),
                     detail::get_or<field::QUAL>(r, std::string_view{}),
                     detail::get_or<field::ID>(r, std::string_view{}),
                     detail::get_or<field::OFFSET>(r, 0u),
                     detail::get_or<field::REF_SEQ>(r, std::string_view{}),
                     detail::get_or<field::REF_ID>(r, std::string_view{}),
                     detail::get_or<field::REF_OFFSET>(r, -1), // 1 is added in format SAM
                     detail::get_or<field::ALIGNMENT>(r, default_align_t{}),
                     detail::get_or<field::MAPQ>(r, 0u),
                     detail::get_or<field::FLAG>(r, 0u),
                     detail::get_or<field::MATE>(r, default_mate_t{}),
                     detail::get_or<field::TAGS>(r, sam_tag_dictionary{}),
                     detail::get_or<field::EVALUE>(r, 0u),
                     detail::get_or<field::BIT_SCORE>(r, 0u));
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
     * Constant. TODO linear in the size of the written alignments?
     *
     * ### Exceptions
     *
     * Basic exception safety.
     *
     * ### Example
     *
     * \snippet test/snippet/io/alignment_file/push_back_tuple.cpp all
     */
    template <typename tuple_t>
    void push_back(tuple_t && t)
    //!\cond
        requires TupleLike<tuple_t>
    //!\endcond
    {
        using default_align_t = std::pair<std::basic_string_view<gapped<char>>, std::basic_string_view<gapped<char>>>;
        using default_mate_t  = std::tuple<std::string_view, uint32_t, uint32_t>;

        // index_of might return npos, but this will be handled well by get_or_ignore (and just return ignore)
        write_record(detail::get_or<selected_field_ids::index_of(field::HEADER_PTR)>(t, nullptr),
                     detail::get_or<selected_field_ids::index_of(field::SEQ)>(t, std::string_view{}),
                     detail::get_or<selected_field_ids::index_of(field::QUAL)>(t, std::string_view{}),
                     detail::get_or<selected_field_ids::index_of(field::ID)>(t, std::string_view{}),
                     detail::get_or<selected_field_ids::index_of(field::OFFSET)>(t, 0u),
                     detail::get_or<selected_field_ids::index_of(field::REF_SEQ)>(t, std::string_view{}),
                     detail::get_or<selected_field_ids::index_of(field::REF_ID)>(t, std::string_view{}),
                     detail::get_or<selected_field_ids::index_of(field::REF_OFFSET)>(t, -1), // 1 is added in format SAM
                     detail::get_or<selected_field_ids::index_of(field::ALIGNMENT)>(t, default_align_t{}),
                     detail::get_or<selected_field_ids::index_of(field::MAPQ)>(t, 0u),
                     detail::get_or<selected_field_ids::index_of(field::FLAG)>(t, 0u),
                     detail::get_or<selected_field_ids::index_of(field::MATE)>(t, default_mate_t{}),
                     detail::get_or<selected_field_ids::index_of(field::TAGS)>(t, sam_tag_dictionary{}),
                     detail::get_or<selected_field_ids::index_of(field::EVALUE)>(t, 0u),
                     detail::get_or<selected_field_ids::index_of(field::BIT_SCORE)>(t, 0u));
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
     * Constant. TODO linear in the size of the written alignments?
     *
     * ### Exceptions
     *
     * Basic exception safety.
     *
     * ### Example
     *
     * \snippet test/snippet/io/alignment_file/emplace_back.cpp all
     */
    template <typename arg_t, typename ... arg_types>
    void emplace_back(arg_t && arg, arg_types && ... args)
    {
        push_back(std::tie(arg, args...));
    }

    /*!\brief            Write a range of records (or tuples) to the file.
     * \tparam rng_t     Type of the range, must satisfy seqan3::OutputRange and have a reference type that
     *                   satisfies seqan3::TupleLike.
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
     * ```cpp
     * alignment_file_output fout{"/tmp/my.sam"};
     *
     * std::vector<std::tuple<dna5_vector, std::string>> range
     * {
     *     { "ACGT"_dna5, "First" },
     *     { "NATA"_dna5, "2nd" },
     *     { "GATA"_dna5, "Third" }
     * }; // a range of "records"
     *
     * fout = range; // will iterate over the records and write them
     * ```
     */
    template <typename rng_t>
    alignment_file_output & operator=(rng_t && range)
    //!\cond
        requires std::ranges::InputRange<rng_t> && TupleLike<reference_t<rng_t>>
    //!\endcond
    {
        for (auto && record : range)
            push_back(std::forward<decltype(record)>(record));
        return *this;
    }

    /*!\brief            Write a range of records (or tuples) to the file.
     * \tparam rng_t     Type of the range, must satisfy seqan3::std::ranges::InputRange and have a reference type that
     *                   satisfies seqan3::TupleLike.
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
     * \snippet test/snippet/io/alignment_file/output.cpp write_range
     *
     * This is especially useful in combination with file-based filters:
     *
     * \todo TODO Implement the snippet once seqan3::alignment_file_in is implemented
     *
     * ```cpp
     * alignment_file_in{"input.sam"} | view::minimum_average_quality_filter(20)
     *                                | view::minimum_alignment_length_filter(50)
     *                                | ranges::view::take(5)
     *                                | alignment_file_output{"output.sam"};
     * ```
     */
    template <typename rng_t>
    friend alignment_file_output & operator|(rng_t && range, alignment_file_output & f)
    //!\cond
        requires std::ranges::InputRange<rng_t> && TupleLike<reference_t<rng_t>>
    //!\endcond
    {
        f = range;
        return f;
    }

    //!\overload
    template <typename rng_t>
    friend alignment_file_output operator|(rng_t && range, alignment_file_output && f)
    //!\cond
        requires std::ranges::InputRange<rng_t> && TupleLike<reference_t<rng_t>>
    //!\endcond
    {
        f = range;
        return std::move(f);
    }
    //!\}

    //!\brief The options are public and its members can be set directly.
    alignment_file_output_options options;

    /*!\cond DEV
     * \brief Expose a reference to the underlying stream object. [public, but not documented as part of the API]
     */
    stream_type & get_stream()
    {
        return stream;
    }
    //!\endcond

    /*!\brief Access the file's header.
     *
     * \details
     *
     * ### Example
     *
     * \snippet test/snippet/io/alignment_file/output.cpp set_header
     *
     * \sa seqan3::alignment_file_header
     */
    alignment_file_header & header()
    {
        if (header_ptr == nullptr)
            header_ptr = std::unique_ptr<alignment_file_header>(new alignment_file_header);

        return *header_ptr;
    }

protected:
    //!\privatesection

    //!\brief The file header object.
    std::unique_ptr<alignment_file_header> header_ptr{nullptr};

    //!\brief Path of the file that the stream operates on.
    std::string file_name;

    //!\brief The stream we are writing to.
    stream_type stream;

    //!\brief Type of the format, an std::variant over the `valid_formats`.
    using format_type = detail::transfer_template_args_onto_t<valid_formats, std::variant>;

    //!\brief The actual std::variant holding a pointer to the detected/selected format.
    format_type format;

    //!\brief Write record to format.
    template <typename ...pack_type>
    void write_record(alignment_file_header const * hdr_ptr, pack_type && ...remainder)
    {
        if (header_ptr == nullptr && hdr_ptr != nullptr)
            header_ptr = std::unique_ptr<alignment_file_header>(new alignment_file_header(*hdr_ptr));

        static_assert((sizeof...(pack_type) == 14), "Wrong parameter list passed to write_record.");

        assert(!format.valueless_by_exception());

        std::visit([&] (auto & f)
        {
            f.write(stream,
                    options,
                    header_ptr,
                    std::forward<pack_type>(remainder)...);

        }, format);
    }

    //!\brief Befriend iterator so it can access the buffers.
    friend iterator;
};

/*!\name Type deduction guides
 * \relates seqan3::alignment_file_output
 * \{
 */
template <Ostream<char>             stream_type,
          AlignmentFileOutputFormat file_format,
          detail::Fields             selected_field_ids>
alignment_file_output(stream_type && _stream, file_format const &, selected_field_ids const &)
    -> alignment_file_output<selected_field_ids,
                         type_list<file_format>,
                         std::remove_reference_t<stream_type>>;
//!\}

} // namespace seqan3

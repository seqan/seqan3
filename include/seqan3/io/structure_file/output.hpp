// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::structure_file_output and corresponding traits classes.
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 */

#pragma once

#include <cassert>
#include <fstream>
#include <optional>
#include <string>
#include <type_traits>
#include <variant>
#include <vector>

#include <range/v3/algorithm/equal.hpp>
#include <range/v3/view/zip.hpp>

#include <seqan3/alphabet/structure/all.hpp>
#include <seqan3/core/type_traits/basic.hpp>
#include <seqan3/core/type_traits/template_inspection.hpp>
#include <seqan3/core/concept/tuple.hpp>
#include <seqan3/io/stream/concept.hpp>
#include <seqan3/io/exception.hpp>
#include <seqan3/std/filesystem>
#include <seqan3/io/record.hpp>
#include <seqan3/io/detail/misc_output.hpp>
#include <seqan3/io/detail/out_file_iterator.hpp>
#include <seqan3/io/detail/record.hpp>
#include <seqan3/io/structure_file/output_format_concept.hpp>
#include <seqan3/io/structure_file/output_options.hpp>
#include <seqan3/io/structure_file/format_vienna.hpp>
#include <seqan3/range/view/convert.hpp>
#include <seqan3/range/view/get.hpp>
#include <seqan3/std/ranges>

namespace seqan3
{

// ----------------------------------------------------------------------------
// structure_file_output
// ----------------------------------------------------------------------------

/*!\brief A class for writing structured sequence files, e.g. Stockholm, Connect, Vienna, ViennaRNA bpp matrix ...
 * \ingroup structure_file
 * \tparam selected_field_ids A seqan3::fields type with the list and order of fields IDs; only relevant if these
 *                            can't be deduced.
 * \tparam valid_formats      A seqan3::type_list of the selectable formats (each must meet
 *                            seqan3::StructureFileOutputFormat).
 * \tparam stream_char_type   The type of the underlying stream device(s); must model seqan3::Char.
 * \details
 *
 * ### Introduction
 *
 * Structured sequence files contain intra-molecular interactions of RNA or protein. Usually, but not necessarily, they
 * contain the nucleotide or amino acid sequences and descriptions as well. Interactions can be represented
 * either as fixed _secondary structure_, where every character is assigned at most one interaction partner
 * (structure of minimum free energy), or an _annotated sequence_, where every character is assigned a set
 * of interaction partners with specific base pair probabilities.
 *
 * The structured sequence file abstraction supports writing ten different fields:
 *
 *   1. seqan3::field::SEQ (sequence)
 *   2. seqan3::field::ID (identifier)
 *   3. seqan3::field::BPP (annotated sequence)
 *   4. seqan3::field::STRUCTURE (secondary structure)
 *   5. seqan3::field::STRUCTURED_SEQ (sequence and structure in one range)
 *   6. seqan3::field::ENERGY (minimum free energy)
 *   7. seqan3::field::REACT (reactivity)
 *   8. seqan3::field::REACT_ERR (reactivity error)
 *   9. seqan3::field::COMMENT (free text)
 *   10. seqan3::field::OFFSET (index of first sequence character)
 *
 * The member functions take any and either of these fields. If the field ID of an argument cannot be deduced, it
 * is assumed to correspond to the field ID of the respective template parameter.
 *
 * ### Construction and specialisation
 *
 * This class comes with two constructors, one for construction from a file name and one for construction from
 * an existing stream and a known format. The first one automatically picks the format based on the extension
 * of the file name. The second can be used if you have a non-file stream, like std::cout or std::ostringstream,
 * that you want to read from and/or if you cannot use file-extension based detection, but know that your output
 * file has a certain format.
 *
 * In most cases the template parameters are deduced completely automatically:
 *
 * \snippet test/snippet/io/structure_file/structure_file_output.cpp temp_param_deduc
 *
 * Writing to std::cout:
 *
 * \snippet test/snippet/io/structure_file/structure_file_output.cpp write_std_out
 *
 * Note that this is not the same as writing `structure_file_output<>` (with angle brackets). In the latter case they are
 * explicitly set to their default values, in the former case
 * [automatic deduction](http://en.cppreference.com/w/cpp/language/class_template_argument_deduction) happens which
 * chooses different parameters depending on the constructor arguments. Prefer deduction over explicit defaults.
 *
 * ### Writing record-wise
 *
 * You can iterate over this file record-wise:
 *
 * \snippet test/snippet/io/structure_file/structure_file_output.cpp iter_by_rec
 *
 * The easiest way to write to a sequence file is to use the push_back() or emplace_back() member functions. These
 * work similarly to how they work on an std::vector. If you pass a tuple to push_back() or give arguments to
 * emplace_back() the seqan3::field ID of the i-th tuple-element/argument is assumed to be the i-th value of
 * selected_field_ids, i.e. by default the first is assumed to be seqan3::field::SEQ, the second seqan3::field::ID
 * and the third one seqan3::field::STRUCTURE. You may give less fields than are selected, if the actual format you are
 * writing to can cope with less
 * (e.g. for Vienna it is sufficient to write seqan3::field::SEQ, seqan3::field::ID and seqan3::field::STRUCTURE,
 * even if selected_field_ids also contains seqan3::field::ENERGY).
 *
 * You may also use the output file's iterator for writing, however, this rarely provides an advantage.
 *
* ### Writing record-wise (custom fields)
 *
 * If you want to pass a combined object for SEQ and STRUCTURE fields to push_back() / emplace_back(), or if you want
 * to change the order of the parameters, you can pass a non-empty fields trait object to the
 * structure_file_output constructor to select the fields that are used for interpreting the arguments.
 *
 * The following snippets demonstrates the usage of such a fields trait object.
 *
 * \snippet test/snippet/io/structure_file/structure_file_output.cpp write_fields
 *
 * A different way of passing custom fields to the file is to pass a seqan3::record – instead of a tuple – to
 * push_back(). The seqan3::record clearly indicates which of its elements has which seqan3::field ID so the file will
 * use that information instead of the template argument. This is especially handy when reading from one file and
 * writing to another, because you don't have to configure the output file to match the input file, it will just work:
 *
 * \snippet test/snippet/io/structure_file/structure_file_output.cpp pass_rec
 *
 * ### Writing record-wise in batches
 *
 * You can write multiple records at once, by assigning to the file:
 *
 * \snippet test/snippet/io/structure_file/structure_file_output.cpp mult_rec
 * ### File I/O pipelines
 *
 * Record-wise writing in batches also works for writing from input files directly to output files, because input
 * files are also input ranges in SeqAn:
 *
 * \snippet test/snippet/io/structure_file/structure_file_output.cpp file_conv
 *
 * This can be combined with file-based views to create I/O pipelines:
 *
 * \snippet test/snippet/io/structure_file/structure_file_output.cpp pipeline
 *
 * ### Column-based writing
 *
 * The record-based interface treats the file as a range of tuples (the records), but in certain situations
 * you might have the data as columns, i.e. a tuple-of-ranges, instead of range-of-tuples.
 *
 * You can use column-based writing in that case, it uses operator=() :
 *
 * \snippet test/snippet/io/structure_file/structure_file_output.cpp data_storage
 * \snippet test/snippet/io/structure_file/structure_file_output.cpp col_based
 *
 * ### Formats
 *
 * Currently, the only implemented format is seqan3::format_vienna. More formats will follow soon.
 */

template <detail::Fields selected_field_ids_ = fields<field::SEQ, field::ID, field::STRUCTURE>,
          detail::TypeListOfStructureFileOutputFormats valid_formats_ = type_list<format_vienna>,
          Char stream_char_type_ = char>
class structure_file_output
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
    //!\brief Character type of the stream(s), usually `char`.
    using stream_char_type      = stream_char_type_;
    //!\}

    //!\brief The subset of seqan3::field IDs that are valid for this file.
    using field_ids = fields<field::SEQ,
                             field::ID,
                             field::BPP,
                             field::STRUCTURE,
                             field::STRUCTURED_SEQ,
                             field::ENERGY,
                             field::REACT,
                             field::REACT_ERR,
                             field::COMMENT,
                             field::OFFSET>;

    static_assert([] () constexpr
                  {
                      for (field f : selected_field_ids::as_array)
                          if (!field_ids::contains(f))
                              return false;
                      return true;
                  }(),
                  "You selected a field that is not valid for structure files, please refer to the documentation "
                  "of structure_file_output::field_ids for the accepted values.");

    static_assert([] () constexpr
                  {
                      return !(selected_field_ids::contains(field::STRUCTURED_SEQ) &&
                               (selected_field_ids::contains(field::SEQ) ||
                               (selected_field_ids::contains(field::STRUCTURE))));
                  }(), "You may not select field::STRUCTURED_SEQ and either of field::SEQ and field::STRUCTURE "
                       "at the same time.");

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
    using iterator          = detail::out_file_iterator<structure_file_output>;
    //!\brief The const iterator type is void, because files are not const-iterable.
    using const_iterator    = void;
    //!\brief The type returned by end().
    using sentinel          = std::ranges::default_sentinel_t;
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    //!\brief Default constructor is explicitly deleted, you need to give a stream or file name.
    structure_file_output()                                          = delete;
    //!\brief Copy construction is explicitly deleted, because you can't have multiple access to the same file.
    structure_file_output(structure_file_output const &)             = delete;
    //!\brief Copy assignment is explicitly deleted, because you can't have multiple access to the same file.
    structure_file_output & operator=(structure_file_output const &) = delete;
    //!\brief Move construction is defaulted.
    structure_file_output(structure_file_output &&)                  = default;
    //!\brief Move assignment is defaulted.
    structure_file_output & operator=(structure_file_output &&)      = default;
    //!\brief Destructor is defaulted.
    ~structure_file_output()                                         = default;

    /*!\brief Construct from filename.
     * \param[in] filename Path to the file you wish to open.
     * \param[in] fields_tag A seqan3::fields tag. [optional]
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
    structure_file_output(std::filesystem::path filename,
                          selected_field_ids const & SEQAN3_DOXYGEN_ONLY(fields_tag) = selected_field_ids{}) :
        primary_stream{new std::ofstream{filename, std::ios_base::out | std::ios::binary}, stream_deleter_default}
    {
        if (!primary_stream->good())
            throw file_open_error{"Could not open file " + filename.string() + " for writing."};

        // possibly add intermediate compression stream
        secondary_stream = detail::make_secondary_ostream(*primary_stream, filename);

        // initialise format handler or throw if format is not found
        detail::set_format(format, filename);
    }

    /*!\brief Construct from an existing stream and with specified format.
     * \tparam file_format The format of the file in the stream, must satisfy
     * seqan3::StructureFileOutputFormat.
     * \param[in,out] stream The stream to write to, must be derived of std::basic_ostream.
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
    template <OStream2 stream_t, StructureFileOutputFormat file_format>
    structure_file_output(stream_t & stream,
                          file_format const & SEQAN3_DOXYGEN_ONLY(format_tag),
                          selected_field_ids const & SEQAN3_DOXYGEN_ONLY(fields_tag) = selected_field_ids{}) :
        primary_stream{&stream, stream_deleter_noop},
        secondary_stream{&stream, stream_deleter_noop},
        format{detail::structure_file_output_format<file_format>{}}
    {
        static_assert(meta::in<valid_formats, file_format>::value,
                      "You selected a format that is not in the valid_formats of this file.");
    }

    //!\overload
    template <OStream2 stream_t, StructureFileOutputFormat file_format>
    structure_file_output(stream_t && stream,
                          file_format const & SEQAN3_DOXYGEN_ONLY(format_tag),
                          selected_field_ids const & SEQAN3_DOXYGEN_ONLY(fields_tag) = selected_field_ids{}) :
        primary_stream{new stream_t{std::move(stream)}, stream_deleter_default},
        secondary_stream{&*primary_stream, stream_deleter_noop},
        format{detail::structure_file_output_format<file_format>{}}
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
     * \snippet test/snippet/io/structure_file/structure_file_output.cpp push_back
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
     * \snippet test/snippet/io/structure_file/structure_file_output.cpp pass_rec
     */
    template <typename record_t>
    void push_back(record_t && r)
        requires TupleLike<record_t> &&
                 requires { requires detail::is_type_specialisation_of_v<remove_cvref_t<record_t>, record>; }
    {
        write_record(detail::get_or_ignore<field::SEQ>(r),
                     detail::get_or_ignore<field::ID>(r),
                     detail::get_or_ignore<field::BPP>(r),
                     detail::get_or_ignore<field::STRUCTURE>(r),
                     detail::get_or_ignore<field::STRUCTURED_SEQ>(r),
                     detail::get_or_ignore<field::ENERGY>(r),
                     detail::get_or_ignore<field::REACT>(r),
                     detail::get_or_ignore<field::REACT_ERR>(r),
                     detail::get_or_ignore<field::COMMENT>(r),
                     detail::get_or_ignore<field::OFFSET>(r));
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
     * \snippet test/snippet/io/structure_file/structure_file_output.cpp push_back_2
     */
    template <typename tuple_t>
    void push_back(tuple_t && t)
        requires TupleLike<tuple_t>
    {
        // index_of might return npos, but this will be handled well by get_or_ignore (and just return ignore)
        write_record(detail::get_or_ignore<selected_field_ids::index_of(field::SEQ)>(t),
                     detail::get_or_ignore<selected_field_ids::index_of(field::ID)>(t),
                     detail::get_or_ignore<selected_field_ids::index_of(field::BPP)>(t),
                     detail::get_or_ignore<selected_field_ids::index_of(field::STRUCTURE)>(t),
                     detail::get_or_ignore<selected_field_ids::index_of(field::STRUCTURED_SEQ)>(t),
                     detail::get_or_ignore<selected_field_ids::index_of(field::ENERGY)>(t),
                     detail::get_or_ignore<selected_field_ids::index_of(field::REACT)>(t),
                     detail::get_or_ignore<selected_field_ids::index_of(field::REACT_ERR)>(t),
                     detail::get_or_ignore<selected_field_ids::index_of(field::COMMENT)>(t),
                     detail::get_or_ignore<selected_field_ids::index_of(field::OFFSET)>(t));
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
     * \snippet test/snippet/io/structure_file/structure_file_output.cpp emplace_back
     */
    template <typename arg_t, typename ... arg_types>
    void emplace_back(arg_t && arg, arg_types && ... args)
    {
        push_back(std::tie(arg, args...));
    }

    /*!\brief            Write a range of records (or tuples) to the file.
     * \tparam rng_t     Type of the range, must satisfy std::ranges::OutputRange and have a reference type that
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
     * \snippet test/snippet/io/structure_file/structure_file_output.cpp equal
     */
    template <std::ranges::InputRange rng_t>
    structure_file_output & operator=(rng_t && range)
        requires TupleLike<reference_t<rng_t>>
    {
        for (auto && record : range)
            push_back(std::forward<decltype(record)>(record));
        return *this;
    }

    /*!\brief            Write a range of records (or tuples) to the file.
     * \tparam rng_t     Type of the range, must satisfy std::ranges::InputRange and have a reference type that
     *                   satisfies seqan3::TupleLike.
     * \param[in] range  The range to write.
     * \param[in] f      The file being written to.
     *
     * \details
     *
     * This operator enables structure_file_output to be at the end of a piping operation. It just calls
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
     * \snippet test/snippet/io/structure_file/structure_file_output.cpp pipe_func
     *
     * This is especially useful in combination with file-based filters:
     *
     * \snippet test/snippet/io/structure_file/structure_file_output.cpp pipeline
     */
    template <std::ranges::InputRange rng_t>
    friend structure_file_output & operator|(rng_t && range, structure_file_output & f)
        requires TupleLike<reference_t<rng_t>>
    {
        f = range;
        return f;
    }

    //!\overload
    template <std::ranges::InputRange rng_t>
    friend structure_file_output operator|(rng_t && range, structure_file_output && f)
        requires TupleLike<reference_t<rng_t>>
    {
        f = range;
        return std::move(f);
    }
    //!\}

    /*!\name Tuple interface
     * \brief Provides functions for field-based ("column"-based) writing.
     * \{
     */
    /*!\brief            Write columns (wrapped in a seqan3::record) to the file.
     * \tparam typelist  Template argument to seqan3::record, each type must be a column (range-of-range).
     * \tparam field_ids Template argument to seqan3::record, the IDs corresponding to the columns.
     * \param[in] r      The record of columns.
     *
     * \details
     *
     * \attention This is not part of the row-based file writing; the seqan3::record does not represent a file record,
     * it is a tuple of the columns (with field information).
     *
     * ### Complexity
     *
     * Linear in the size of the columns.
     *
     * ### Exceptions
     *
     * Basic exception safety.
     *
     * ### Example
     *
     * \snippet test/snippet/io/structure_file/structure_file_output.cpp col_based
     */
    template <typename typelist, typename field_ids>
    structure_file_output & operator=(record<typelist, field_ids> const & r)
    {
        write_columns(detail::range_wrap_ignore(detail::get_or_ignore<field::SEQ>(r)),
                      detail::range_wrap_ignore(detail::get_or_ignore<field::ID>(r)),
                      detail::range_wrap_ignore(detail::get_or_ignore<field::BPP>(r)),
                      detail::range_wrap_ignore(detail::get_or_ignore<field::STRUCTURE>(r)),
                      detail::range_wrap_ignore(detail::get_or_ignore<field::STRUCTURED_SEQ>(r)),
                      detail::range_wrap_ignore(detail::get_or_ignore<field::ENERGY>(r)),
                      detail::range_wrap_ignore(detail::get_or_ignore<field::REACT>(r)),
                      detail::range_wrap_ignore(detail::get_or_ignore<field::REACT_ERR>(r)),
                      detail::range_wrap_ignore(detail::get_or_ignore<field::COMMENT>(r)),
                      detail::range_wrap_ignore(detail::get_or_ignore<field::OFFSET>(r)));
        return *this;
    }

    /*!\brief            Write columns (wrapped in a std::tuple) to the file.
     * \tparam arg_types The column types, each type must be a range-of-range.
     * \param[in] t      The tuple of columns.
     *
     * \details
     *
     * The columns are assumed to correspond to the field IDs given in selected_field_ids, however passing less
     * is accepted if the format does not require all of them.
     *
     * ### Complexity
     *
     * Linear in the size of the columns.
     *
     * ### Exceptions
     *
     * Basic exception safety.
     *
     * ### Example
     *
     * \snippet test/snippet/io/structure_file/structure_file_output.cpp col_based
     *
     */
    template <typename ... arg_types>
    structure_file_output & operator=(std::tuple<arg_types...> const & t)
    {
        // index_of might return npos, but this will be handled well by get_or_ignore (and just return ignore)
        write_columns(
           detail::range_wrap_ignore(detail::get_or_ignore<selected_field_ids::index_of(field::SEQ)>(t)),
           detail::range_wrap_ignore(detail::get_or_ignore<selected_field_ids::index_of(field::ID)>(t)),
           detail::range_wrap_ignore(detail::get_or_ignore<selected_field_ids::index_of(field::BPP)>(t)),
           detail::range_wrap_ignore(detail::get_or_ignore<selected_field_ids::index_of(field::STRUCTURE)>(t)),
           detail::range_wrap_ignore(detail::get_or_ignore<selected_field_ids::index_of(field::STRUCTURED_SEQ)>(t)),
           detail::range_wrap_ignore(detail::get_or_ignore<selected_field_ids::index_of(field::ENERGY)>(t)),
           detail::range_wrap_ignore(detail::get_or_ignore<selected_field_ids::index_of(field::REACT)>(t)),
           detail::range_wrap_ignore(detail::get_or_ignore<selected_field_ids::index_of(field::REACT_ERR)>(t)),
           detail::range_wrap_ignore(detail::get_or_ignore<selected_field_ids::index_of(field::COMMENT)>(t)),
           detail::range_wrap_ignore(detail::get_or_ignore<selected_field_ids::index_of(field::OFFSET)>(t)));
        return *this;
    }
    //!\}

    //!\brief The options are public and its members can be set directly.
    structure_file_output_options options;

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
    using format_type = typename detail::variant_from_tags<valid_formats, detail::structure_file_output_format>::type;
    //!\brief The actual std::variant holding a pointer to the detected/selected format.
    format_type format;
    //!\}

    //!\brief Write record to format.
    template <typename seq_type,
              typename id_type,
              typename bpp_type,
              typename structure_type,
              typename structured_seq_type,
              typename energy_type,
              typename react_type,
              typename comment_type,
              typename offset_type>
    void write_record(seq_type && seq,
                      id_type && id,
                      bpp_type && bpp,
                      structure_type && structure,
                      structured_seq_type && structured_seq,
                      energy_type && energy,
                      react_type && react,
                      react_type && react_error,
                      comment_type && comment,
                      offset_type && offset)
    {
        static_assert(detail::decays_to_ignore_v<structured_seq_type> ||
                      (detail::decays_to_ignore_v<seq_type> && detail::decays_to_ignore_v<structure_type>),
                      "You may not select field::STRUCTURED_SEQ and either of field::SEQ and field::STRUCTURE "
                      "at the same time.");

        assert(!format.valueless_by_exception());
        std::visit([&] (auto & f)
        {
            if constexpr (!detail::decays_to_ignore_v<structured_seq_type>)
            {
                f.write(*secondary_stream,
                        options,
                        structured_seq | view::get<0>,
                        id,
                        bpp,
                        structured_seq | view::get<1>,
                        energy,
                        react,
                        react_error,
                        comment,
                        offset);
            }
            else
            {
                f.write(*secondary_stream,
                        options,
                        seq,
                        id,
                        bpp,
                        structure,
                        energy,
                        react,
                        react_error,
                        comment,
                        offset);
            }
        }, format);
    }

    //!\brief Write columns to file format, only tag-dispatch once.
    template <std::ranges::InputRange seq_type,
              std::ranges::InputRange id_type,
              std::ranges::InputRange bpp_type,
              std::ranges::InputRange structure_type,
              std::ranges::InputRange structured_seq_type,
              std::ranges::InputRange energy_type,
              std::ranges::InputRange react_type,
              std::ranges::InputRange comment_type,
              std::ranges::InputRange offset_type>
    void write_columns(seq_type && seq,
                       id_type && id,
                       bpp_type && bpp,
                       structure_type && structure,
                       structured_seq_type && structured_seq,
                       energy_type && energy,
                       react_type && react,
                       react_type && react_error,
                       comment_type && comment,
                       offset_type && offset)
    {
        static_assert(!(detail::decays_to_ignore_v<reference_t<seq_type>> &&
                        detail::decays_to_ignore_v<reference_t<id_type>> &&
                        detail::decays_to_ignore_v<reference_t<bpp_type>> &&
                        detail::decays_to_ignore_v<reference_t<structure_type>> &&
                        detail::decays_to_ignore_v<reference_t<structured_seq_type>> &&
                        detail::decays_to_ignore_v<reference_t<energy_type>> &&
                        detail::decays_to_ignore_v<reference_t<react_type>> &&
                        detail::decays_to_ignore_v<reference_t<comment_type>> &&
                        detail::decays_to_ignore_v<reference_t<offset_type>>),
                      "At least one of the columns must not be set to std::ignore.");

        static_assert(detail::decays_to_ignore_v<reference_t<structured_seq_type>> ||
                      (detail::decays_to_ignore_v<reference_t<seq_type>> &&
                       detail::decays_to_ignore_v<reference_t<structure_type>>),
                      "You may not select field::STRUCTURED_SEQ and either of field::SEQ and field::STRUCTURE "
                      "at the same time.");

        assert(!format.valueless_by_exception());
        std::visit([&] (auto & f)
        {
            if constexpr (!detail::decays_to_ignore_v<reference_t<structured_seq_type>>)
            {
                auto zipped = std::view::zip(structured_seq, id, bpp, energy, react, react_error, comment, offset);

                for (auto && v : zipped)
                {
                    f.write(*secondary_stream,
                            options,
                            std::get<0>(v) | view::get<0>, // seq
                            std::get<1>(v),  // id
                            std::get<2>(v),  // bpp
                            std::get<0>(v) | view::get<1>, // structure
                            std::get<3>(v),  // energy
                            std::get<4>(v),  // react
                            std::get<5>(v),  // react_error
                            std::get<6>(v),  // comment
                            std::get<7>(v)); // offset
                }
            }
            else
            {
                auto zipped = std::view::zip(seq, id, bpp, structure, energy, react, react_error, comment, offset);

                for (auto && v : zipped)
                {
                    f.write(*secondary_stream, options, std::get<0>(v), std::get<1>(v), std::get<2>(v), std::get<3>(v),
                            std::get<4>(v), std::get<5>(v), std::get<6>(v), std::get<7>(v), std::get<8>(v));
                }
            }
        }, format);
    }

    //!\brief Befriend iterator so it can access the buffers.
    friend iterator;
};

/*!\name Type deduction guides
 * \relates seqan3::structure_file_output
 * \{
 */

//!\brief Deduction of the selected fields, the file format and the stream type.
template <OStream2                  stream_t,
          StructureFileOutputFormat file_format,
          detail::Fields            selected_field_ids>
structure_file_output(stream_t &&, file_format const &, selected_field_ids const &)
    -> structure_file_output<selected_field_ids,
                             type_list<file_format>,
                             typename std::remove_reference_t<stream_t>::char_type>;

//!\overload
template <OStream2                  stream_t,
          StructureFileOutputFormat file_format,
          detail::Fields            selected_field_ids>
structure_file_output(stream_t &, file_format const &, selected_field_ids const &)
    -> structure_file_output<selected_field_ids,
                             type_list<file_format>,
                             typename std::remove_reference_t<stream_t>::char_type>;
//!\}

} // namespace seqan3

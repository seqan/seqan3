// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin 
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik 
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::structure_file_out and corresponding traits classes.
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
#include <seqan3/core/metafunction/basic.hpp>
#include <seqan3/core/metafunction/template_inspection.hpp>
#include <seqan3/core/concept/tuple.hpp>
#include <seqan3/io/stream/concept.hpp>
#include <seqan3/io/exception.hpp>
#include <seqan3/io/filesystem.hpp>
#include <seqan3/io/record.hpp>
#include <seqan3/io/detail/out_file_iterator.hpp>
#include <seqan3/io/detail/record.hpp>
#include <seqan3/io/structure_file/output_format_concept.hpp>
#include <seqan3/io/structure_file/output_options.hpp>
#include <seqan3/io/structure_file/format_vienna.hpp>
#include <seqan3/range/view/convert.hpp>
#include <seqan3/std/ranges>

namespace seqan3
{

// ----------------------------------------------------------------------------
// structure_file_out
// ----------------------------------------------------------------------------

/*!\brief A class for writing structured sequence files, e.g. Stockholm, Connect, Vienna, ViennaRNA bpp matrix ...
 * \ingroup structure_file
 * \tparam selected_field_ids A seqan3::fields type with the list and order of fields IDs; only relevant if these
 *                            can't be deduced.
 * \tparam valid_formats      A seqan3::type_list of the selectable formats (each must meet
 *                            seqan3::structure_file_output_format_concept).
 * \tparam stream_type        The type of the stream, must satisfy seqan3::ostream_concept.
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
 * ```cpp
 * structure_file_out fout{"/tmp/my.dbn"}; // Vienna format detected, std::ofstream opened for file
 * ```
 *
 * Writing to std::cout:
 * ```cpp
 * structure_file_out fout{std::move(std::cout), structure_file_format_vienna{}};
 * //               ^ no need to specify the template arguments
 *
 * fout.emplace_back("example_id", "ACGTN"_dna5);
 * ```
 *
 * Note that this is not the same as writing `structure_file_out<>` (with angle brackets). In the latter case they are
 * explicitly set to their default values, in the former case
 * [automatic deduction](http://en.cppreference.com/w/cpp/language/class_template_argument_deduction) happens which
 * chooses different parameters depending on the constructor arguments. For opening from file, `structure_file_out<>`
 * would have also worked, but for opening from stream it would not have.
 *
 * ### Writing record-wise
 *
 * You can iterate over this file record-wise:
 *
 * ```cpp
 * structure_file_out fout{"/tmp/my.dbn"};
 *
 * for // ...
 * {
 *     std::string id;
 *     rna5_vector seq;
 *     std::vector<wuss51> structure;
 *
 *     // ...
 *
 *     fout.emplace_back(seq, id, structure);        // as individual variables
 *     // or:
 *     fout.push_back(std::tie(seq, id, structure)); // as a tuple
 * }
 * ```
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
 * structure_file_out constructor to select the fields that are used for interpreting the arguments.
 *
 * The following snippets demonstrates the usage of such a fields trait object.
 *
 * ```cpp
 * structure_file_out fout{"/tmp/my.dbn", fields<field::ID, field::STRUCTURED_SEQ>{}};
 *
 * for // ...
 * {
 *     std::string id;
 *     std::vector<structured_rna<rna5, wuss51>> structured_seq; // vector of combined data structure
 *
 *     // ...
 *
 *     fout.emplace_back(id, structured_seq); // note also that the order the argumets is now different, because
 *     // or:                                    you specified that ID should be first in the fields template argument
 *     fout.push_back(std::tie(id, structured_seq));
 * }
 * ```
 *
 * A different way of passing custom fields to the file is to pass a seqan3::record – instead of a tuple – to
 * push_back(). The seqan3::record clearly indicates which of its elements has which seqan3::field ID so the file will
 * use that information instead of the template argument. This is especially handy when reading from one file and
 * writing to another, because you don't have to configure the output file to match the input file, it will just work:
 *
 * ```cpp
 * structure_file_in  fin{"input.dbn", fields<field::ID, field::STRUCTURED_SEQ>{}};
 * structure_file_out fout{"output.dbn"}; // doesn't have to match the configuration
 *
 * for (auto & r : fin)
 * {
 *     if // r fulfills some criterium
 *     {
 *         fout.push_back(r);
 *     }
 * }
 * ```
 *
 * ### Writing record-wise in batches
 *
 * You can write multiple records at once, by assigning to the file:
 *
 * ```cpp
 * structure_file_out fout{"/tmp/my.dbn"};
 *
 * std::vector<std::tuple<rna5_vector, std::string, std::vector<wuss51>>> range
 * {
 *     { "ACGT"_rna5, "First", "...."_wuss51 },
 *     { "NATA"_rna5, "2nd",   "...."_wuss51 },
 *     { "GATA"_rna5, "Third", "...."_wuss51 }
 * }; // a range of "records"
 *
 * fout = range; // will iterate over the records and write them
 * ```
 * ### File I/O pipelines
 *
 * Record-wise writing in batches also works for writing from input files directly to output files, because input
 * files are also input ranges in SeqAn:
 *
 * ```cpp
 * // file format conversion in one line:
 * structure_file_out fout{"output.dbn"} = structure_file_in{"input.dbn"};
 *
 * // or in pipe notation:
 * structure_file_in{"input.dbn"} | structure_file_out{"output.dbn"};
 * ```
 *
 * This can be combined with file-based views to create I/O pipelines:
 *
 * ```cpp
 * structure_file_in{"input.dbn"} | view::minimum_sequence_length_filter(50)
 *                                | ranges::view::take(5)
 *                                | structure_file_out{"output.dbn"};
 * ```
 *
 * ### Column-based writing
 *
 * The record-based interface treats the file as a range of tuples (the records), but in certain situations
 * you might have the data as columns, i.e. a tuple-of-ranges, instead of range-of-tuples.
 *
 * You can use column-based writing in that case, it uses operator=() :
 *
 * ```cpp
 *
 * struct data_storage_t
 * {
 *     concatenated_sequences<rna5_vector>         sequences;
 *     concatenated_sequences<std::string>         ids;
 *     concatenated_sequences<std::vector<wuss51>> structures;
 * };
 *
 * data_storage_t data_storage; // a global or globally used variable in your program
 *
 * // ... in your file writing function:
 *
 * structure_file_out fout{"/tmp/my.dbn"};
 *
 * fout = std::tie(data_storage.sequences, data_storage.ids, data_storage.structures);
 * ```
 *
 * ### Formats
 *
 * Currently, the only implemented format is seqan3::structure_file_format_vienna. More formats will follow soon.
 */

template <detail::fields_concept selected_field_ids_ = fields<field::SEQ, field::ID, field::STRUCTURE>,
          detail::type_list_of_structure_file_output_formats_concept valid_formats_
              = type_list<structure_file_format_vienna>,
          ostream_concept<char> stream_type_ = std::ofstream>
class structure_file_out
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
                  "of structure_file_out::field_ids for the accepted values.");

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
    using value_type        = void;
    using reference         = void;
    using const_reference   = void;
    using size_type         = void;
    //!\brief A signed integer type, usually std::ptrdiff_t.
    using difference_type   = std::ptrdiff_t;
    //!\brief The iterator type of this view (an output iterator).
    using iterator          = detail::out_file_iterator<structure_file_out>;
    //!\brief The const iterator type is void, because files are not const-iterable.
    using const_iterator    = void;
    //!\brief The type returned by end().
    using sentinel          = std::ranges::default_sentinel;
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    //!\brief Default constructor is explicitly deleted, you need to give a stream or file name.
    structure_file_out() = delete;
    //!\brief Copy construction is explicitly deleted, because you can't have multiple access to the same file.
    structure_file_out(structure_file_out const &) = delete;
    //!\brief Copy assignment is explicitly deleted, because you can't have multiple access to the same file.
    structure_file_out & operator=(structure_file_out const &) = delete;
    //!\brief Move construction is defaulted.
    structure_file_out(structure_file_out &&) = default;
    //!\brief Move assignment is defaulted.
    structure_file_out & operator=(structure_file_out &&) = default;
    //!\brief Destructor is defaulted.
    ~structure_file_out() = default;

    /*!\brief Construct from filename.
     * \param[in] _file_name Path to the file you wish to open.
     * \param[in] fields_tag A seqan3::fields tag. [optional]
     *
     * \details
     *
     * In addition to the file name, you may specify a custom seqan3::fields type which may be easier than
     * defining all the template parameters.
     */
    structure_file_out(filesystem::path const & _file_name,
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
     * \tparam file_format The format of the file in the stream, must satisfy
     * seqan3::structure_file_output_format_concept.
     * \param[in] _stream  The stream to operate on (this must be std::move'd in!).
     * \param[in] format_tag The file format tag.
     * \param[in] fields_tag A seqan3::fields tag. [optional]
     */
    template <structure_file_output_format_concept file_format>
    structure_file_out(stream_type             && _stream,
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
     * ```cpp
     * structure_file_out fout{"/tmp/my.dbn"};
     *
     * auto it = fout.begin();
     *
     * for // ...
     * {
     *     std::string id;
     *     rna5_vector seq;
     *     std::vector<wuss51> structure;
     *
     *     // ...
     *
     *     // assign to iterator
     *     *it = std::tie(seq, id, structure);
     *     // is the same as:
     *     fout.push_back(std::tie(seq, id, structure));
     * }
     * ```
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
     * ```cpp
     * structure_file_out fout{"/tmp/my.dbn"};
     *
     * auto it = fout.begin();
     *
     * for // ...
     * {
     *     std::string id;
     *     rna5_vector seq;
     *     std::vector<wuss51> structure;
     *
     *     // ...
     *
     *     fout.push_back(std::tie(seq, id, structure));
     * }
     * ```
     */
    template <typename tuple_t>
    void push_back(tuple_t && t)
        requires tuple_like_concept<tuple_t>
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
     * ```cpp
     * structure_file_out fout{"/tmp/my.dbn"};
     *
     * auto it = fout.begin();
     *
     * for // ...
     * {
     *     std::string id;
     *     rna5_vector seq;
     *     std::vector<wuss51> structure;
     *
     *     // ...
     *
     *     fout.emplace_back(seq, id, structure);
     * }
     * ```
     */
    template <typename arg_t, typename ... arg_types>
    void emplace_back(arg_t && arg, arg_types && ... args)
    {
        push_back(std::tie(arg, args...));
    }

    /*!\brief            Write a range of records (or tuples) to the file.
     * \tparam rng_t     Type of the range, must satisfy std::ranges::OutputRange and have a reference type that
     *                   satisfies seqan3::tuple_like_concept.
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
     * structure_file_out fout{"/tmp/my.dbn"};
     *
     * std::vector<std::tuple<rna5_vector, std::string, std::vector<wuss51>>> range
     * {
     *     { "ACGT"_rna5, "First", "...."_wuss51 },
     *     { "NATA"_rna5, "2nd",   "...."_wuss51 },
     *     { "GATA"_rna5, "Third", "...."_wuss51 }
     * }; // a range of "records"
     *
     * fout = range; // will iterate over the records and write them
     * ```
     */
    template <std::ranges::InputRange rng_t>
    structure_file_out & operator=(rng_t && range)
        requires tuple_like_concept<reference_t<rng_t>>
    {
        for (auto && record : range)
            push_back(std::forward<decltype(record)>(record));
        return *this;
    }

    /*!\brief            Write a range of records (or tuples) to the file.
     * \tparam rng_t     Type of the range, must satisfy std::ranges::InputRange and have a reference type that
     *                   satisfies seqan3::tuple_like_concept.
     * \param[in] range  The range to write.
     * \param[in] f      The file being written to.
     *
     * \details
     *
     * This operator enables structure_file_out to be at the end of a piping operation. It just calls
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
     * ```cpp
     * structure_file_out fout{"/tmp/my.dbn"};
     *
     * std::vector<std::tuple<rna5_vector, std::string, std::vector<wuss51>>> range
     * {
     *     { "ACGT"_rna5, "First", "...."_wuss51 },
     *     { "NATA"_rna5, "2nd",   "...."_wuss51 },
     *     { "GATA"_rna5, "Third", "...."_wuss51 }
     * }; // a range of "records"
     *
     * range | fout;
     * // the same as:
     * fout = range;
     * ```
     *
     * This is especially useful in combination with file-based filters:
     *
     * ```cpp
     * structure_file_in{"input.dbn"} | view::minimum_sequence_length_filter(50)
     *                                | ranges::view::take(5)
     *                                | structure_file_out{"output.dbn"};
     * ```
     */
    template <std::ranges::InputRange rng_t>
    friend structure_file_out & operator|(rng_t && range, structure_file_out & f)
        requires tuple_like_concept<reference_t<rng_t>>
    {
        f = range;
        return f;
    }

    //!\overload
    template <std::ranges::InputRange rng_t>
    friend structure_file_out operator|(rng_t && range, structure_file_out && f)
        requires tuple_like_concept<reference_t<rng_t>>
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
     * ```cpp
     *
     * struct data_storage_t
     * {
     *     concatenated_sequences<rna5_vector>         sequences;
     *     concatenated_sequences<std::string>         ids;
     *     concatenated_sequences<std::vector<wuss51>> structures;
     * };
     *
     * data_storage_t data_storage; // a global or globally used variable in your program
     *
     * // ... in your file writing function:
     *
     * structure_file_out fout{"/tmp/my.dbn"};
     *
     * fout = std::tie(data_storage.sequences, data_storage.ids, data_storage.structures);
     * ```
     */
    template <typename typelist, typename field_ids>
    structure_file_out & operator=(record<typelist, field_ids> const & r)
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
     * ```cpp
     *
     * struct data_storage_t
     * {
     *     concatenated_sequences<rna5_vector>         sequences;
     *     concatenated_sequences<std::string>         ids;
     *     concatenated_sequences<std::vector<wuss51>> structures;
     * };
     *
     * data_storage_t data_storage; // a global or globally used variable in your program
     *
     * // ... in your file writing function:
     *
     * structure_file_out fout{"/tmp/my.dbn"};
     *
     * fout = std::tie(data_storage.sequences, data_storage.ids, data_storage.structures);
     * ```
     */
    template <typename ... arg_types>
    structure_file_out & operator=(std::tuple<arg_types...> const & t)
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
     * \brief Expose a reference to the underlying stream object. [public, but not documented as part of the API]
     */
    stream_type & get_stream()
    {
        return stream;
    }
    //!\endcond
protected:
    //!\privatesection
    //!\brief Path of the file that the stream operates on.
    std::string file_name;

    //!\brief The stream we are writing to.
    stream_type stream;

    //!\brief Type of the format, an std::variant over the `valid_formats`.
    using format_type = detail::transfer_template_args_onto_t<valid_formats, std::variant>;
    //!\brief The actual std::variant holding a pointer to the detected/selected format.
    format_type format;

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
                f.write(stream,
                        options,
                        structured_seq | view::convert<typename structured_seq_type::sequence_alphabet_type>,
                        id,
                        bpp,
                        structured_seq | view::convert<typename structured_seq_type::structure_alphabet_type>,
                        energy,
                        react,
                        react_error,
                        comment,
                        offset);
            }
            else
            {
                f.write(stream,
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
                auto zipped = ranges::view::zip(structured_seq, id, bpp, energy, react, react_error, comment, offset);

                for (auto && v : zipped)
                {
                    f.write(stream,
                            options,
                            std::get<0>(v) | view::convert
                                             <typename reference_t<structured_seq_type>::sequence_alphabet_type>,
                            std::get<1>(v),  // id
                            std::get<2>(v),  // bpp
                            std::get<0>(v) | view::convert
                                             <typename reference_t<structured_seq_type>::structure_alphabet_type>,
                            std::get<3>(v),  // energy
                            std::get<4>(v),  // react
                            std::get<5>(v),  // react_error
                            std::get<6>(v),  // comment
                            std::get<7>(v)); // offset
                }
            }
            else
            {
                auto zipped = ranges::view::zip(seq, id, bpp, structure, energy, react, react_error, comment, offset);

                for (auto && v : zipped)
                {
                    f.write(stream, options, std::get<0>(v), std::get<1>(v), std::get<2>(v), std::get<3>(v),
                            std::get<4>(v), std::get<5>(v), std::get<6>(v), std::get<7>(v), std::get<8>(v));
                }
            }
        }, format);
    }

    //!\brief Befriend iterator so it can access the buffers.
    friend iterator;
};

/*!\name Type deduction guides
 * \relates seqan3::structure_file_out
 * \{
 */
template <ostream_concept<char>                stream_type,
          structure_file_output_format_concept file_format,
          detail::fields_concept               selected_field_ids>
structure_file_out(stream_type && _stream, file_format const &, selected_field_ids const &)
    -> structure_file_out<selected_field_ids,
                         type_list<file_format>,
                         std::remove_reference_t<stream_type>>;
//!\}

} // namespace seqan3

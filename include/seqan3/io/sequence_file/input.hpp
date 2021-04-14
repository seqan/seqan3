// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::sequence_file_input and corresponding traits classes.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <cassert>
#include <seqan3/std/filesystem>
#include <fstream>
#include <string>
#include <variant>
#include <vector>

#include <seqan3/alphabet/adaptation/char.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/nucleotide/dna15.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/alphabet/quality/qualified.hpp>
#include <seqan3/core/detail/pack_algorithm.hpp>
#include <seqan3/io/detail/in_file_iterator.hpp>
#include <seqan3/io/detail/misc_input.hpp>
#include <seqan3/io/detail/record.hpp>
#include <seqan3/io/exception.hpp>
#include <seqan3/io/sam_file/format_sam.hpp>
#include <seqan3/io/sequence_file/format_embl.hpp>
#include <seqan3/io/sequence_file/format_fasta.hpp>
#include <seqan3/io/sequence_file/format_fastq.hpp>
#include <seqan3/io/sequence_file/format_genbank.hpp>
#include <seqan3/io/sequence_file/input_format_concept.hpp>
#include <seqan3/io/sequence_file/record.hpp>
#include <seqan3/io/stream/concept.hpp>
#include <seqan3/utility/type_list/traits.hpp>

namespace seqan3
{

// ----------------------------------------------------------------------------
// sequence_file_input_traits
// ----------------------------------------------------------------------------

/*!\interface seqan3::sequence_file_input_traits <>
 * \brief The requirements a traits_type for seqan3::sequence_file_input must meet.
 * \ingroup io_sequence_file
 */
/*!\name Requirements for seqan3::sequence_file_input_traits
 * \brief You can expect these **member types** of all types that satisfy seqan3::sequence_file_input_traits.
 * \memberof seqan3::sequence_file_input_traits
 *
 * \details
 *
 * \{
 */
/*!\typedef using sequence_alphabet
 * \brief Alphabet of the characters for the seqan3::field::seq; must satisfy seqan3::alphabet.
 */
/*!\typedef using sequence_legal_alphabet
 * \brief Intermediate alphabet for seqan3::field::seq; must satisfy seqan3::alphabet and be convertible to
 * `sequence_alphabet`.
 *
 * \details
 *
 * This alphabet can be a superset of `sequence_alphabet` to allow conversion of some characters
 * without producing an error, e.g. if this is set to seqan3::dna15 and `sequence_alphabet` is set to seqan3::dna5,
 * 'M' will be an accepted character and automatically converted to 'N', while 'Z' will still be an illegal
 * character and produce an error.
 */
/*!\typedef using sequence_container
 * \brief Type template of the seqan3::field::seq, a container template over `sequence_alphabet`;
 * must satisfy seqan3::sequence_container.
 */
/*!\typedef using id_alphabet
 * \brief Alphabet of the characters for the seqan3::field::id; must satisfy seqan3::alphabet.
 */
/*!\typedef using id_container
 * \brief Type template of the seqan3::field::id, a container template over `id_alphabet`;
 * must satisfy seqan3::sequence_container.
 */
/*!\typedef using quality_alphabet
 * \brief Alphabet of the characters for the seqan3::field::qual; must satisfy seqan3::writable_quality_alphabet.
 */
/*!\typedef using quality_container
 * \brief Type template of the seqan3::field::qual, a container template over `quality_alphabet`;
 * must satisfy seqan3::sequence_container.
 */
//!\}
//!\cond
template <typename t>
SEQAN3_CONCEPT sequence_file_input_traits = requires (t v)
{
    requires writable_alphabet<typename t::sequence_alphabet>;
    requires writable_alphabet<typename t::sequence_legal_alphabet>;
    requires explicitly_convertible_to<typename t::sequence_legal_alphabet, typename t::sequence_alphabet>;
    requires sequence_container<typename t::template sequence_container<typename t::sequence_alphabet>>;

    requires writable_alphabet<typename t::id_alphabet>;
    requires sequence_container<typename t::template id_container<typename t::id_alphabet>>;

    requires writable_quality_alphabet<typename t::quality_alphabet>;
    requires sequence_container<typename t::template quality_container<typename t::quality_alphabet>>;
};
//!\endcond

// ----------------------------------------------------------------------------
// sequence_file_input_default_traits
// ----------------------------------------------------------------------------

/*!\brief The default traits for seqan3::sequence_file_input
 * \implements sequence_file_input_traits
 * \ingroup io_sequence_file
 *
 * \details
 *
 * If you wish to change a single or a few types from the default, just inherit from this class and
 * "overwrite" the respective type definitions.
 *
 * This example will make the file read into a smaller alphabet and a compressed container:
 *
 * \include test/snippet/io/sequence_file/sequence_file_input_trait_overwrite.cpp
 */
struct sequence_file_input_default_traits_dna
{
    /*!\name Member types
     * \brief Definitions to satisfy seqan3::sequence_file_input_traits.
     * \{
     */

    //!\brief The sequence alphabet is seqan3::dna5.
    using sequence_alphabet                 = dna5;

    //!\brief The legal sequence alphabet for parsing is seqan3::dna15.
    using sequence_legal_alphabet           = dna15;

    //!\brief The type of a DNA sequence is std::vector.
    template <typename _sequence_alphabet>
    using sequence_container                = std::vector<_sequence_alphabet>;

    //!\brief The alphabet for an identifier string is char.
    using id_alphabet                       = char;

    //!\brief The string type for an identifier is std::basic_string.
    template <typename _id_alphabet>
    using id_container                      = std::basic_string<_id_alphabet>;

    //!\brief The alphabet for a quality annotation is seqan3::phred42.
    using quality_alphabet                  = phred42;

    //!\brief The string type for a quality annotation is std::vector.
    template <typename _quality_alphabet>
    using quality_container                 = std::vector<_quality_alphabet>;

    //!\}
};

//!\brief A traits type that specifies input as amino acids.
//!\ingroup io_sequence_file
struct sequence_file_input_default_traits_aa : sequence_file_input_default_traits_dna
{
    /*!\name Member types
     * \brief Definitions to satisfy seqan3::sequence_file_input_traits.
     * \{
     */

    //!\brief The sequence alphabet is seqan3::aa27.
    using sequence_alphabet = aa27;

    //!\brief The legal sequence alphabet for parsing is seqan3::aa27.
    using sequence_legal_alphabet = aa27;
    //!\}
};

// ----------------------------------------------------------------------------
// sequence_file_input
// ----------------------------------------------------------------------------

/*!\brief A class for reading sequence files, e.g. FASTA, FASTQ ...
 * \ingroup io_sequence_file
 * \tparam traits_type          An auxiliary type that defines certain member types and constants, must satisfy
 *                              seqan3::sequence_file_input_traits.
 * \tparam selected_field_ids   A seqan3::fields type with the list and order of desired record entries; all fields
 *                              must be in seqan3::sequence_file_input::field_ids.
 * \tparam valid_formats        A seqan3::type_list of the selectable formats (each must meet
 *                              seqan3::sequence_file_input_format).
 *
 * \details
 *
 * ### Introduction
 *
 * Sequence files are the most generic and common biological files. Well-known formats include
 * FastA and FastQ, but some may also be interested in treating SAM or BAM files as sequence
 * files, discarding the alignment.
 *
 * The Sequence file abstraction supports reading three different fields:
 *
 *   1. seqan3::field::seq
 *   2. seqan3::field::id
 *   3. seqan3::field::qual
 *
 * The first three fields are retrieved by default (and in that order). The last field may be selected to have
 * sequence and qualities directly stored in a more memory-efficient combined container. If you select the last
 * field you may not select seqan3::field::seq or seqan3::field::qual.
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
 * \snippet test/snippet/io/sequence_file/sequence_file_input_template_deduction.cpp main
 * Reading from an std::istringstream:
 * \include test/snippet/io/sequence_file/sequence_file_input_istringstream.cpp
 *
 * Note that this is not the same as writing `sequence_file_input<>` (with angle brackets). In the latter case they are
 * explicitly set to their default values, in the former case
 * [automatic deduction](https://en.cppreference.com/w/cpp/language/class_template_argument_deduction) happens which
 * chooses different parameters depending on the constructor arguments. For opening from file, `sequence_file_input<>`
 * would have also worked, but for opening from stream it would not have.
 *
 * In some cases, you do need to specify the arguments, e.g. if you want to read amino acids:
 *
 * \include test/snippet/io/sequence_file/sequence_file_input_aminoacid.cpp
 *
 * You can define your own traits type to further customise the types used by and returned by this class, see
 * seqan3::sequence_file_input_default_traits_dna for more details. As mentioned above, specifying at least one
 * template parameter yourself means that you loose automatic deduction so if you want to read amino acids **and**
 * want to read from a string stream you need to give all types yourself:
 *
 * \include test/snippet/io/sequence_file/sequence_file_input_template_specification.cpp
 *
 * ### Reading record-wise
 *
 * You can iterate over this file record-wise:
 *
 * \include test/snippet/io/sequence_file/sequence_file_input_record_iter.cpp
 *
 * In the above example, `record` has the type \ref record_type which is seqan3::sequence_record.
 *
 * *Note:* It is important to write `auto &` and not just `auto`, otherwise you will copy the record on every iteration.
 * Since the buffer gets "refilled" on every iteration, you can also move the data out of the record if you want
 * to store it somewhere without copying:
 *
 * \include test/snippet/io/sequence_file/sequence_file_input_auto_ref.cpp
 *
 * ### Reading record-wise (decomposed records)
 *
 * Instead of using member accessor on the record, you can also use
 * [structured bindings](https://en.cppreference.com/w/cpp/language/structured_binding)
 * to decompose the record into its elements:
 *
 * \include test/snippet/io/sequence_file/sequence_file_input_decomposed.cpp
 *
 * In this case you immediately get the two elements of the tuple: `sequence` of \ref sequence_type and `id` of
 * \ref id_type. **But beware: with structured bindings you do need to get the order of elements correctly!**
 *
 * ### Reading record-wise (custom fields)
 *
 * If you want to skip specific fields from the record you can pass a non-empty fields trait object to the
 * sequence_file_input constructor to select the fields that should be read from the input. For example to choose a
 * combined field for SEQ and QUAL (see above). Or to never actually read the QUAL, if you don't need it.
 * The following snippets demonstrate the usage of such a fields trait object.
 *
 * \include test/snippet/io/sequence_file/sequence_file_input_custom_fields.cpp
 *
 * When reading a file, all fields not present in the file (but requested implicitly or via the `selected_field_ids`
 * parameter) are ignored.
 *
 * ### Views on files
 *
 * Since SeqAn files are ranges, you can also create views over files. A useful example is to filter the records
 * based on certain criteria, e.g. minimum length of the sequence field:
 *
 * \include test/snippet/io/sequence_file/sequence_file_input_file_view.cpp
 *
 * ### End of file
 *
 * You can check whether a file is at end by comparing begin() and end() (if they are the same, the file is at end).
 *
 * ### Formats
 *
 * We currently support reading the following formats:
 *   * seqan3::format_fasta
 *   * seqan3::format_fastq
 *   * seqan3::format_embl
 *   * seqan3::format_genbank
 *   * seqan3::format_sam
 */

template <
    sequence_file_input_traits traits_type_ = sequence_file_input_default_traits_dna,
    detail::fields_specialisation selected_field_ids_ = fields<field::seq, field::id, field::qual>,
    detail::type_list_of_sequence_file_input_formats valid_formats_ = type_list<format_embl,
                                                                                format_fasta,
                                                                                format_fastq,
                                                                                format_genbank,
                                                                                format_sam>>
class sequence_file_input
{
public:
    /*!\name Template arguments
     * \brief Exposed as member types for public access.
     * \{
     */
    //!\brief A traits type that defines aliases and template for storage of the fields.
    using traits_type           = traits_type_;
    //!\brief A seqan3::fields list with the fields selected for the record.
    using selected_field_ids    = selected_field_ids_;
    //!\brief A seqan3::type_list with the possible formats.
    using valid_formats         = valid_formats_;
    //!\brief Character type of the stream(s).
    using stream_char_type      = char;
    //!\}

    /*!\brief The subset of seqan3::field IDs that are valid for this file; order corresponds to the types in
     * \ref field_types.
     */
    using field_ids            = fields<field::seq, field::id, field::qual>;

    static_assert([] () constexpr
                  {
                      for (field f : selected_field_ids::as_array)
                          if (!field_ids::contains(f))
                              return false;
                      return true;
                  }(),
                  "You selected a field that is not valid for sequence files, please refer to the documentation "
                  "of sequence_file_input::field_ids for the accepted values.");

    /*!\name Field types and record type
     * \brief These types are relevant for record/row-based reading; they may be manipulated via the \ref traits_type
     * to achieve different storage behaviour.
     * \{
     */
    //!\brief The type of field::seq (std::vector <seqan3::dna5> by default).
    using sequence_type         = typename traits_type::template sequence_container<
                                    typename traits_type::sequence_alphabet>;
    //!\brief The type of field::id (std::string by defaul).
    using id_type               = typename traits_type::template id_container<
                                    typename traits_type::id_alphabet>;
    //!\brief The type of field::qual (std::vector <seqan3::phred42> by default).
    using quality_type          = typename traits_type::template quality_container<
                                    typename traits_type::quality_alphabet>;
    //!\brief The previously defined types aggregated in a seqan3::type_list.
    using field_types           = type_list<sequence_type, id_type, quality_type>;

    //!\brief The type of the record, a specialisation of seqan3::record; acts as a tuple of the selected field types.
    using record_type           = sequence_record<detail::select_types_with_ids_t<field_types,
                                                                                  field_ids,
                                                                                  selected_field_ids>,
                                                  selected_field_ids>;
    //!\}

    /*!\name Range associated types
     * \brief The types necessary to facilitate the behaviour of an input range (used in record-wise reading).
     * \{
     */
    //!\brief The value_type is the \ref record_type.
    using value_type        = record_type;
    //!\brief The reference type.
    using reference         = record_type &;
    //!\brief The const_reference type is void, because files are not const-iterable.
    using const_reference   = void;
    //!\brief An unsigned integer type, usually std::size_t.
    using size_type         = size_t;
    //!\brief A signed integer type, usually std::ptrdiff_t.
    using difference_type   = std::make_signed_t<size_t>;
    //!\brief The iterator type of this view (an input iterator).
    using iterator          = detail::in_file_iterator<sequence_file_input>;
    //!\brief The const iterator type is void, because files are not const-iterable.
    using const_iterator    = void;
    //!\brief The type returned by end().
    using sentinel          = std::default_sentinel_t;
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    //!\brief Default constructor is explicitly deleted, you need to give a stream or file name.
    sequence_file_input() = delete;
    //!\brief Copy construction is explicitly deleted, because you can't have multiple access to the same file.
    sequence_file_input(sequence_file_input const &) = delete;
    //!\brief Copy assignment is explicitly deleted, because you can't have multiple access to the same file.
    sequence_file_input & operator=(sequence_file_input const &) = delete;
    //!\brief Move construction is defaulted.
    sequence_file_input(sequence_file_input &&) = default;
    //!\brief Move assignment is defaulted.
    sequence_file_input & operator=(sequence_file_input &&) = default;
    //!\brief Destructor is defaulted.
    ~sequence_file_input() = default;

    /*!\brief Construct from filename.
     * \param[in] filename      Path to the file you wish to open.
     * \param[in] fields_tag    A seqan3::fields tag. [optional]
     * \throws seqan3::file_open_error If the file could not be opened, e.g. non-existant, non-readable, unknown format.
     *
     * \details
     *
     * In addition to the file name, you may specify a custom seqan3::fields type which may be easier than
     * defining all the template parameters.
     *
     * ### Decompression
     *
     * This constructor transparently applies a decompression stream on top of the file stream in case
     * the file is detected as being compressed.
     * See the section on \link io_compression compression and decompression \endlink for more information.
     */
    sequence_file_input(std::filesystem::path filename,
                        selected_field_ids const & SEQAN3_DOXYGEN_ONLY(fields_tag) = selected_field_ids{}) :
        primary_stream{new std::ifstream{}, stream_deleter_default}
    {
        primary_stream->rdbuf()->pubsetbuf(stream_buffer.data(), stream_buffer.size());
        static_cast<std::basic_ifstream<char> *>(primary_stream.get())->open(filename,
                                                                             std::ios_base::in | std::ios::binary);

        if (!primary_stream->good())
            throw file_open_error{"Could not open file " + filename.string() + " for reading."};

        // possibly add intermediate compression stream
        secondary_stream = detail::make_secondary_istream(*primary_stream, filename);

        // initialise format handler or throw if format is not found
        using format_variant_t = typename detail::variant_from_tags<valid_formats,
                                                                    detail::sequence_file_input_format_exposer>::type;
        format_variant_t format_variant{};
        detail::set_format(format_variant, filename);

        std::visit([&] (auto && selected_format)
        {
            using format_t = std::remove_cvref_t<decltype(selected_format)>;
            format = std::make_unique<selected_sequence_format<format_t>>();
        }, format_variant);
    }
    /* NOTE(h-2): Curiously we do not need a user-defined deduction guide for the above constructor.
     * A combination of default template parameters and auto-deduction guides works as expected,
     * independent of whether the second/optional parameter is specified or not, i.e. it is possible
     * to auto-deduct and overwrite a single template parameter out of the four if the optional parameter
     * is specified and use the default otherwise.
     */

    /*!\brief Construct from an existing stream and with specified format.
     * \tparam file_format   The format of the file in the stream, must satisfy seqan3::sequence_file_input_format.
     * \param[in] stream     The stream to operate on; must be derived of std::basic_istream.
     * \param[in] format_tag The file format tag.
     * \param[in] fields_tag A seqan3::fields tag. [optional]
     *
     * \details
     *
     * ### Decompression
     *
     * This constructor transparently applies a decompression stream on top of the stream in case
     * it is detected as being compressed.
     * See the section on \link io_compression compression and decompression \endlink for more information.
     */
    template <input_stream stream_t,
              sequence_file_input_format file_format>
    //!\cond
        requires std::same_as<typename std::remove_reference_t<stream_t>::char_type, stream_char_type>
    //!\endcond
    sequence_file_input(stream_t                 & stream,
                        file_format        const & SEQAN3_DOXYGEN_ONLY(format_tag),
                        selected_field_ids const & SEQAN3_DOXYGEN_ONLY(fields_tag) = selected_field_ids{}) :
        primary_stream{&stream, stream_deleter_noop},
        format{std::make_unique<selected_sequence_format<file_format>>()}
    {
        static_assert(list_traits::contains<file_format, valid_formats>,
                      "You selected a format that is not in the valid_formats of this file.");

        // possibly add intermediate compression stream
        secondary_stream = detail::make_secondary_istream(*primary_stream);
    }

    //!\overload
    template <input_stream stream_t,
              sequence_file_input_format file_format>
    //!\cond
        requires std::same_as<typename std::remove_reference_t<stream_t>::char_type, stream_char_type>
    //!\endcond
    sequence_file_input(stream_t                && stream,
                        file_format        const & SEQAN3_DOXYGEN_ONLY(format_tag),
                        selected_field_ids const & SEQAN3_DOXYGEN_ONLY(fields_tag) = selected_field_ids{}) :
        primary_stream{new stream_t{std::move(stream)}, stream_deleter_default},
        format{std::make_unique<selected_sequence_format<file_format>>()}
    {
        static_assert(list_traits::contains<file_format, valid_formats>,
                      "You selected a format that is not in the valid_formats of this file.");

        // possibly add intermediate compression stream
        secondary_stream = detail::make_secondary_istream(*primary_stream);
    }
    //!\}

    /*!\name Range interface
     * \brief Provides functions for record based reading of the file.
     * \{
     */
    /*!\brief Returns an iterator to current position in the file.
     * \returns An iterator pointing to the current position in the file.
     * \throws seqan3::format_error
     *
     * Equals end() if the file is at end.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * Throws seqan3::format_error if the first record could not be read into the buffer.
     */
    iterator begin()
    {
        // buffer first record
        if (!first_record_was_read)
        {
            read_next_record();
            first_record_was_read = true;
        }

        return {*this};
    }

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

    /*!\brief Return the record we are currently at in the file.
     * \returns A reference to the currently buffered record.
     *
     * This function returns a reference to the currently buffered record, it is identical to dereferencing begin(),
     * but begin also always points to the current record on single pass input ranges:
     *
     * \include test/snippet/io/sequence_file/sequence_file_input_return_record.cpp
     *
     * It most situations using the iterator interface or a range-based for-loop are preferable to using front(),
     * because you can only move to the next record via the iterator.
     *
     * In any case, don't forget the reference! If you want to save the data from the record elsewhere, use move:
     *
     * \include test/snippet/io/sequence_file/sequence_file_input_record_move.cpp
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    reference front() noexcept
    {
        return *begin();
    }
    //!\}

    //!\brief The input file options type.
    using sequence_file_input_options_type = sequence_file_input_options<typename traits_type::sequence_legal_alphabet>;
    //!\brief The options are public and its members can be set directly.
    sequence_file_input_options_type options;

protected:
    //!\privatesection
    /*!\name Data buffers
     * \{
     */
    //!\brief Buffer for a single record.
    record_type record_buffer;
    //!\brief A larger (compared to stl default) stream buffer to use when reading from a file.
    std::vector<char> stream_buffer{std::vector<char>(1'000'000)};
    //!\brief Buffer for the previous record position.
    std::streampos position_buffer{};
    //!\}

    /*!\name Stream / file access
     * \{
     */
    //!\brief The type of the internal stream pointers. Allows dynamically setting ownership management.
    using stream_ptr_t = std::unique_ptr<std::basic_istream<stream_char_type>,
                                         std::function<void(std::basic_istream<stream_char_type>*)>>;
    //!\brief Stream deleter that does nothing (no ownership assumed).
    static void stream_deleter_noop(std::basic_istream<stream_char_type> *) {}
    //!\brief Stream deleter with default behaviour (ownership assumed).
    static void stream_deleter_default(std::basic_istream<stream_char_type> * ptr) { delete ptr; }

    //!\brief The primary stream is the user provided stream or the file stream if constructed from filename.
    stream_ptr_t primary_stream{nullptr, stream_deleter_noop};
    //!\brief The secondary stream is a compression layer on the primary or just points to the primary (no compression).
    stream_ptr_t secondary_stream{nullptr, stream_deleter_noop};

    //!\brief Tracks whether the very first record is buffered when calling begin().
    bool first_record_was_read{false};
    //!\brief File is at position 1 behind the last record.
    bool at_end{false};
    //!\}

private:
    //!\brief Tell the format to move to the next record and update the buffer.
    void read_next_record()
    {
        // clear the record
        record_buffer.clear();

        // at end if we could not read further
        if ((std::istreambuf_iterator<stream_char_type>{*secondary_stream} ==
             std::istreambuf_iterator<stream_char_type>{}))
        {
            at_end = true;
            return;
        }

        format->read_sequence_record(*secondary_stream, record_buffer, position_buffer, options);
    }

    /*!\brief An abstract base class to store the selected input format.
     *
     * \details
     *
     * This abstract base class is used to store the user given input format as a type-erased base class in a
     * std::unique_ptr. There is exactly one derived type that implements this abstract base class called
     * seqan3::sequence_file_input::selected_sequence_format which is a private subclass of
     * seqan3::sequence_file_input. It is not exposed to the user and allows to hide the implementation detail of
     * storing a specfic format instance which is first known at runtime.
     */
    struct sequence_format_base
    {
        /*!\name Constructors, destructor and assignment
         * \{
         */
        sequence_format_base() = default; //!< Default.
        sequence_format_base(sequence_format_base const &) = default; //!< Default.
        sequence_format_base(sequence_format_base &&) = default; //!< Default.
        sequence_format_base & operator=(sequence_format_base const &) = default; //!< Default.
        sequence_format_base & operator=(sequence_format_base &&) = default; //!< Default.
        virtual ~sequence_format_base() = default; //!< Virtual default.
        //!\}

        /*!\brief Reads the next format specific record from the given istream.
         *
         * \param[in, out] instream The input stream to extract the next record from.
         * \param[in, out] record_buffer The record buffer to fill.
         * \param[in, out] position_buffer The buffer to store the position of the current record.
         * \param[in] options User specific format options set from outside.
         *
         * \details
         *
         * Invokes the actual read sequence record function for the selected format and fills the record accordingly.
         */
        virtual void read_sequence_record(std::istream & instream,
                                          record_type & record_buffer,
                                          std::streampos & position_buffer,
                                          sequence_file_input_options_type const & options) = 0;
    };

    /*!\brief The specific selected format to read the records from.
     *
     * \tparam format_t The user specific format type to store.
     *
     * \details
     *
     * This class implements the format specific read operation based on the instantiated format type.
     * The specfic format type is selected at runtime and stored through a pointer to the abstract
     * seqan3::sequence_file_input::sequence_format_base class. A virtual function call then ensures that the
     * specific read record function of the selected format is invoked.
     */
    template <typename format_t>
    struct selected_sequence_format final : public sequence_format_base
    {
        /*!\name Constructors, destructor and assignment
         * \{
         */
        selected_sequence_format() = default; //!< Default.
        selected_sequence_format(selected_sequence_format const &) = default; //!< Default.
        selected_sequence_format(selected_sequence_format &&) = default; //!< Default.
        selected_sequence_format & operator=(selected_sequence_format const &) = default; //!< Default.
        selected_sequence_format & operator=(selected_sequence_format &&) = default; //!< Default.
        ~selected_sequence_format() = default; //!< Default.
        //!\}

        //!\copydoc sequence_format_base::read_sequence_record
        void read_sequence_record(std::istream & instream,
                                  record_type & record_buffer,
                                  std::streampos & position_buffer,
                                  sequence_file_input_options_type const & options) override
        {
            // read new record
            {
                _format.read_sequence_record(instream,
                                             options,
                                             position_buffer,
                                             detail::get_or_ignore<field::seq>(record_buffer),
                                             detail::get_or_ignore<field::id>(record_buffer),
                                             detail::get_or_ignore<field::qual>(record_buffer));
            }
        };

        //!\brief The selected format stored as a format exposer object.
        detail::sequence_file_input_format_exposer<format_t> _format{};
    };

    //!\brief An instance of the detected/selected format.
    std::unique_ptr<sequence_format_base> format{};

    //!\brief Befriend iterator so it can access the buffers.
    friend iterator;
};

/*!\name Type deduction guides
 * \relates seqan3::sequence_file_input
 * \{
 */

//!\brief Deduces the sequence input file type from the stream and the format.
template <input_stream stream_type,
          sequence_file_input_format file_format>
sequence_file_input(stream_type & stream,
                    file_format const &)
    -> sequence_file_input<typename sequence_file_input<>::traits_type,         // actually use the default
                           typename sequence_file_input<>::selected_field_ids,  // default field ids.
                           type_list<file_format>>;

//!\overload
template <input_stream stream_type,
          sequence_file_input_format file_format>
sequence_file_input(stream_type && stream,
                    file_format const &)
    -> sequence_file_input<typename sequence_file_input<>::traits_type,         // actually use the default
                           typename sequence_file_input<>::selected_field_ids,  // default field ids.
                           type_list<file_format>>;

//!\brief Deduces the sequence input file type from the stream, the format and the field ids.
template <input_stream stream_type,
          sequence_file_input_format file_format,
          detail::fields_specialisation selected_field_ids>
sequence_file_input(stream_type && stream,
                    file_format const &,
                    selected_field_ids const &)
    -> sequence_file_input<typename sequence_file_input<>::traits_type,       // actually use the default
                           selected_field_ids,
                           type_list<file_format>>;

//!\overload
template <input_stream stream_type,
          sequence_file_input_format file_format,
          detail::fields_specialisation selected_field_ids>
sequence_file_input(stream_type & stream,
                    file_format const &,
                    selected_field_ids const &)
    -> sequence_file_input<typename sequence_file_input<>::traits_type,       // actually use the default
                           selected_field_ids,
                           type_list<file_format>>;
//!\}

} // namespace seqan3

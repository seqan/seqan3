// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::sequence_file_input and corresponding traits classes.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <cassert>
#include <fstream>
#include <string>
#include <variant>
#include <vector>

#include <range/v3/algorithm/equal.hpp>

#include <seqan3/alphabet/adaptation/char.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/alphabet/quality/qualified.hpp>
#include <seqan3/core/type_traits/basic.hpp>
#include <seqan3/io/stream/concept.hpp>
#include <seqan3/io/exception.hpp>
#include <seqan3/std/filesystem>
#include <seqan3/io/record.hpp>
#include <seqan3/io/detail/in_file_iterator.hpp>
#include <seqan3/io/detail/misc_input.hpp>
#include <seqan3/io/detail/record.hpp>
#include <seqan3/io/sequence_file/input_format_concept.hpp>
#include <seqan3/io/sequence_file/format_embl.hpp>
#include <seqan3/io/sequence_file/format_fasta.hpp>
#include <seqan3/io/sequence_file/format_fastq.hpp>
#include <seqan3/io/sequence_file/format_genbank.hpp>
#include <seqan3/io/sequence_file/format_sam.hpp>
#include <seqan3/range/container/concatenated_sequences.hpp>

namespace seqan3
{

// ----------------------------------------------------------------------------
// SequenceFileInputTraits
// ----------------------------------------------------------------------------

/*!\interface seqan3::SequenceFileInputTraits <>
 * \brief The requirements a traits_type for seqan3::sequence_file_input must meet.
 * \ingroup sequence
 */
/*!\name Requirements for seqan3::SequenceFileInputTraits
 * \brief You can expect these **member types** of all types that satisfy seqan3::SequenceFileInputTraits.
 * \memberof seqan3::SequenceFileInputTraits
 *
 * \details
 *
 * Note that the alphabet type of the seqan3::field::SEQ_QUAL cannot be specified directly, it is always
 * seqan3::qualified<sequence_alphabet, quality_alphabet> and the container type templates for
 * the field are those of seqan3::field::SEQ.
 *
 * \{
 */
/*!\typedef using sequence_alphabet
 * \brief Alphabet of the characters for the seqan3::field::SEQ; must satisfy seqan3::Alphabet.
 */
/*!\typedef using sequence_legal_alphabet
 * \brief Intermediate alphabet for seqan3::field::SEQ; must satisfy seqan3::Alphabet and be convertible to
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
 * \brief Type template of the seqan3::field::SEQ, a container template over `sequence_alphabet`;
 * must satisfy seqan3::SequenceContainer.
 */
/*!\typedef using sequence_container_container
 * \brief Type template of a column of seqan3::field::SEQ, a container template that can hold multiple
 * `sequence_container`; must satisfy seqan3::SequenceContainer.
 */
/*!\typedef using id_alphabet
 * \brief Alphabet of the characters for the seqan3::field::ID; must satisfy seqan3::Alphabet.
 */
/*!\typedef using id_container
 * \brief Type template of the seqan3::field::ID, a container template over `id_alphabet`;
 * must satisfy seqan3::SequenceContainer.
 */
/*!\typedef using id_container_container
 * \brief Type template of a column of seqan3::field::ID, a container template that can hold multiple
 * `id_container`; must satisfy seqan3::SequenceContainer.
 */
/*!\typedef using quality_alphabet
 * \brief Alphabet of the characters for the seqan3::field::QUAL; must satisfy seqan3::WritableQualityAlphabet.
 */
/*!\typedef using quality_container
 * \brief Type template of the seqan3::field::QUAL, a container template over `quality_alphabet`;
 * must satisfy seqan3::SequenceContainer.
 */
/*!\typedef using quality_container_container
 * \brief Type template of a column of seqan3::field::QUAL, a container template that can hold multiple
 * `quality_container`; must satisfy seqan3::SequenceContainer.
 */
//!\}
//!\cond
template <typename t>
SEQAN3_CONCEPT SequenceFileInputTraits = requires (t v)
{
    requires WritableAlphabet<typename t::sequence_alphabet>;
    requires WritableAlphabet<typename t::sequence_legal_alphabet>;
    requires ExplicitlyConvertibleTo<typename t::sequence_legal_alphabet, typename t::sequence_alphabet>;
    requires SequenceContainer<typename t::template sequence_container<typename t::sequence_alphabet>>;
    requires SequenceContainer<typename t::template sequence_container_container<
        typename t::template sequence_container<typename t::sequence_alphabet>>>;

    requires WritableAlphabet<typename t::id_alphabet>;
    requires SequenceContainer<typename t::template id_container<typename t::id_alphabet>>;
    requires SequenceContainer<typename t::template id_container_container<typename t::template id_container<
        typename t::id_alphabet>>>;

    requires WritableQualityAlphabet<typename t::quality_alphabet>;
    requires SequenceContainer<typename t::template quality_container<typename t::quality_alphabet>>;
    requires SequenceContainer<typename t::template quality_container_container<
        typename t::template quality_container<typename t::quality_alphabet>>>;
};
//!\endcond

// ----------------------------------------------------------------------------
// sequence_file_input_default_traits
// ----------------------------------------------------------------------------

/*!\brief The default traits for seqan3::sequence_file_input
 * \implements SequenceFileInputTraits
 * \ingroup sequence
 *
 * \details
 *
 * If you wish to change a single or a few types from the default, just inherit from this class and
 * "overwrite" the respective type definitions.
 *
 * This example will make the file read into a smaller alphabet and a compressed container:
 *
 * \snippet test/snippet/io/sequence_file/sequence_file_input.cpp trait_overwrite
 */
struct sequence_file_input_default_traits_dna
{
    /*!\name Member types
     * \brief Definitions to satisfy seqan3::SequenceFileInputTraits.
     * \{
     */

    //!\brief The sequence alphabet is seqan3::dna5.
    using sequence_alphabet                 = dna5;

    //!\brief The legal sequence alphabet for parsing is seqan3::dna15.
    using sequence_legal_alphabet           = dna15;

    //!\brief The type of a DNA sequence is std::vector.
    template <typename _sequence_alphabet>
    using sequence_container                = std::vector<_sequence_alphabet>;

    //!\brief The container for sequences is seqan3::concatenated_sequences.
    template <typename _sequence_container>
    using sequence_container_container      = concatenated_sequences<_sequence_container>;

    //!\brief The alphabet for an identifier string is char.
    using id_alphabet                       = char;

    //!\brief The string type for an identifier is std::basic_string.
    template <typename _id_alphabet>
    using id_container                      = std::basic_string<_id_alphabet>;

    //!\brief The container for identifier strings is seqan3::concatenated_sequences.
    template <typename _id_container>
    using id_container_container            = concatenated_sequences<_id_container>;

    //!\brief The alphabet for a quality annotation is seqan3::phred42.
    using quality_alphabet                  = phred42;

    //!\brief The string type for a quality annotation is std::vector.
    template <typename _quality_alphabet>
    using quality_container                 = std::vector<_quality_alphabet>;

    //!\brief The container for quality annotation strings is seqan3::concatenated_sequences.
    template <typename _quality_container>
    using quality_container_container       = concatenated_sequences<_quality_container>;
    //!\}
};

//!\brief A traits type that specifies input as amino acids.
//!\ingroup sequence
struct sequence_file_input_default_traits_aa : sequence_file_input_default_traits_dna
{
    /*!\name Member types
     * \brief Definitions to satisfy seqan3::SequenceFileInputTraits.
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
 * \ingroup sequence
 * \tparam traits_type          An auxiliary type that defines certain member types and constants, must satisfy
 * seqan3::SequenceFileInputTraits.
 * \tparam selected_field_ids   A seqan3::fields type with the list and order of desired record entries; all fields
 * must be in seqan3::sequence_file_input::field_ids.
 * \tparam valid_formats        A seqan3::type_list of the selectable formats (each must meet
 * seqan3::SequenceFileInputFormat).
 * \tparam stream_char_type     The type of the underlying stream device(s); must model seqan3::Char.
 * \details
 *
 * ### Introduction
 *
 * Sequence files are the most generic and common biological files. Well-known formats include
 * FastA and FastQ, but some may also be interested in treating SAM or BAM files as sequence
 * files, discarding the alignment.
 *
 * The Sequence file abstraction supports reading four different fields:
 *
 *   1. seqan3::field::SEQ
 *   2. seqan3::field::ID
 *   3. seqan3::field::QUAL
 *   4. seqan3::field::SEQ_QUAL (sequence and qualities in one range)
 *
 * The first three fields are retrieved by default (and in that order). The last field may be selected to have
 * sequence and qualities directly stored in a more memory-efficient combined container. If you select the last
 * field you may not select seqan3::field::SEQ or seqan3::field::QUAL.
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
 * \snippet test/snippet/io/sequence_file/sequence_file_input.cpp template_deduction
 * Reading from an std::istringstream:
 * \snippet test/snippet/io/sequence_file/sequence_file_input.cpp istringstream
 *
 * Note that this is not the same as writing `sequence_file_input<>` (with angle brackets). In the latter case they are
 * explicitly set to their default values, in the former case
 * [automatic deduction](http://en.cppreference.com/w/cpp/language/class_template_argument_deduction) happens which
 * chooses different parameters depending on the constructor arguments. For opening from file, `sequence_file_input<>`
 * would have also worked, but for opening from stream it would not have.
 *
 * In some cases, you do need to specify the arguments, e.g. if you want to read amino acids:
 *
 * \snippet test/snippet/io/sequence_file/sequence_file_input.cpp aminoacid
 *
 * You can define your own traits type to further customise the types used by and returned by this class, see
 * seqan3::sequence_file_default_traits_dna for more details. As mentioned above, specifying at least one
 * template parameter yourself means that you loose automatic deduction so if you want to read amino acids **and**
 * want to read from a string stream you need to give all types yourself:
 *
 * \snippet test/snippet/io/sequence_file/sequence_file_input.cpp template_specification
 *
 * ### Reading record-wise
 *
 * You can iterate over this file record-wise:
 *
 * \snippet test/snippet/io/sequence_file/sequence_file_input.cpp record_iter
 *
 * In the above example, rec has the type \ref record_type which is a specialisation of seqan3::record and behaves
 * like an std::tuple (that's why we can access it via get). Instead of using the seqan3::field based interface on
 * the record, you could also use `std::get<0>` or even `std::get<dna4_vector>` to retrieve the sequence, but it is
 * not recommended, because it is more error-prone.
 *
 * *Note:* It is important to write `auto &` and not just `auto`, otherwise you will copy the record on every iteration.
 * Since the buffer gets "refilled" on every iteration, you can also move the data out of the record if you want
 * to store it somewhere without copying:
 *
 * \snippet test/snippet/io/sequence_file/sequence_file_input.cpp auto_ref
 *
 * ### Reading record-wise (decomposed records)
 *
 * Instead of using `get` on the record, you can also use
 * [structured bindings](http://en.cppreference.com/w/cpp/language/structured_binding)
 * to decompose the record into its elements:
 *
 * \snippet test/snippet/io/sequence_file/sequence_file_input.cpp decomposed
 *
 * In this case you immediately get the two elements of the tuple: `seq` of \ref sequence_type and `id` of
 * \ref id_type. **But beware: with structured bindings you do need to get the order of elements correctly!**
 *
 * ### Reading record-wise (custom fields)
 *
 * If you want to skip specific fields from the record you can pass a non-empty fields trait object to the
 * sequence_file_input constructor to select the fields that should be read from the input. For example to choose a
 * combined field for SEQ and QUAL (see above). Or to never actually read the QUAL, if you don't need it.
 * The following snippets demonstrate the usage of such a fields trait object.
 *
 * \snippet test/snippet/io/sequence_file/sequence_file_input.cpp custom_fields
 *
 * When reading a file, all fields not present in the file (but requested implicitly or via the `selected_field_ids`
 * parameter) are ignored.
 *
 * ### Views on files
 *
 * Since SeqAn files are ranges, you can also create views over files. A useful example is to filter the records
 * based on certain criteria, e.g. minimum length of the sequence field:
 *
 * \snippet test/snippet/io/sequence_file/sequence_file_input.cpp file_view
 *
 * ### End of file
 *
 * You can check whether a file is at end by comparing begin() and end() (if they are the same, the file is at end).
 *
 * ### Column-based reading
 *
 * The record-based interface treats the file as a range of tuples (the records), but in certain situations it
 * is desirable to read the file by field, i.e. column wise (tuple-of-ranges, instead of range-of-tuples).
 *
 * This interface is less flexible, but can save you copy operations in certain scenarios, given that
 * you have sufficient memory to load the entire file at once:
 *
 * \snippet test/snippet/io/sequence_file/sequence_file_input.cpp data_storage
 * \snippet test/snippet/io/sequence_file/sequence_file_input.cpp col_read
 *
 * Note that for this to make sense, your storage data types need to be identical to the corresponding column types
 * of the file. If you require different column types you can specify you own traits, see
 * seqan3::SequenceFileInputTraits.
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
    SequenceFileInputTraits                    traits_type_        = sequence_file_input_default_traits_dna,
    detail::Fields                             selected_field_ids_ = fields<field::SEQ,
                                                                            field::ID,
                                                                            field::QUAL>,
    detail::TypeListOfSequenceFileInputFormats valid_formats_      = type_list<format_embl,
                                                                               format_fasta,
                                                                               format_fastq,
                                                                               format_genbank,
                                                                               format_sam>,
    Char                                       stream_char_type_   = char>
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
    //!\brief Character type of the stream(s), usually `char`.
    using stream_char_type      = stream_char_type_;
    //!\}

    /*!\brief The subset of seqan3::field IDs that are valid for this file; order corresponds to the types in
     * \ref field_types.
     */
    using field_ids            = fields<field::SEQ, field::ID, field::QUAL, field::SEQ_QUAL>;

    static_assert([] () constexpr
                  {
                      for (field f : selected_field_ids::as_array)
                          if (!field_ids::contains(f))
                              return false;
                      return true;
                  }(),
                  "You selected a field that is not valid for sequence files, please refer to the documentation "
                  "of sequence_file_input::field_ids for the accepted values.");

    static_assert([] () constexpr
                  {
                      return !(selected_field_ids::contains(field::SEQ_QUAL) &&
                               (selected_field_ids::contains(field::SEQ) ||
                               (selected_field_ids::contains(field::QUAL))));
                  }(),
                  "You may not select field::SEQ_QUAL and either of field::SEQ and field::QUAL at the same time.");

    /*!\name Field types and record type
     * \brief These types are relevant for record/row-based reading; they may be manipulated via the \ref traits_type
     * to achieve different storage behaviour.
     * \{
     */
    //!\brief The type of field::SEQ (std::vector <seqan3::dna5> by default).
    using sequence_type         = typename traits_type::template sequence_container<
                                    typename traits_type::sequence_alphabet>;
    //!\brief The type of field::ID (std::string by defaul).
    using id_type               = typename traits_type::template id_container<
                                    typename traits_type::id_alphabet>;
    //!\brief The type of field::QUAL (std::vector <seqan3::phred42> by default).
    using quality_type          = typename traits_type::template quality_container<
                                    typename traits_type::quality_alphabet>;
    //!\brief The type of field::SEQ_QUAL (std::vector <seqan3::dna5q> by default).
    using sequence_quality_type = typename traits_type::
                                    template sequence_container<qualified<typename traits_type::sequence_alphabet,
                                                                          typename traits_type::quality_alphabet>>;

    //!\brief The previously defined types aggregated in a seqan3::type_list.
    using field_types           = type_list<sequence_type, id_type, quality_type, sequence_quality_type>;

    //!\brief The type of the record, a specialisation of seqan3::record; acts as a tuple of the selected field types.
    using record_type           = record<detail::select_types_with_ids_t<field_types, field_ids, selected_field_ids>,
                                         selected_field_ids>;
    //!\}

    /*!\name Field column types and tuple type
     * \brief These types are relevant for field/column-wise reading; they may be manipulated via the \ref traits_type
     * to achieve different storage behaviour.
     * \{
     */
    //!\brief Column type of field::SEQ (seqan3::concatenated_sequences<sequence_type> by default).
    using sequence_column_type          = typename traits_type::template sequence_container_container<sequence_type>;
    //!\brief Column type of field::ID (seqan3::concatenated_sequences<id_type> by default).
    using id_column_type                = typename traits_type::template id_container_container<id_type>;
    //!\brief Column type of field::QUAL (seqan3::concatenated_sequences<quality_type> by default).
    using quality_column_type           = typename traits_type::template quality_container_container<quality_type>;
    //!\brief Column type of field::SEQ_QUAL (seqan3::concatenated_sequences<sequence_quality_type> by default).
    using sequence_quality_column_type  = typename traits_type::template sequence_container_container<sequence_quality_type>;
    //!\brief The previously defined types aggregated in a seqan3::type_list.
    using field_column_types            = type_list<sequence_column_type,
                                                    id_column_type,
                                                    quality_column_type,
                                                    sequence_quality_column_type>;
    //!\brief The type emulated by the file when read column-wise.
    using file_as_tuple_type            = record<detail::select_types_with_ids_t<field_column_types,
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
    using sentinel          = std::ranges::default_sentinel_t;
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
        primary_stream{new std::ifstream{filename, std::ios_base::in | std::ios::binary}, stream_deleter_default}
    {
        if (!primary_stream->good())
            throw file_open_error{"Could not open file " + filename.string() + " for reading."};

        // possibly add intermediate compression stream
        secondary_stream = detail::make_secondary_istream(*primary_stream, filename);

        // initialise format handler or throw if format is not found
        detail::set_format(format, filename);

        // buffer first record
        read_next_record();
    }
    /* NOTE(h-2): Curiously we do not need a user-defined deduction guide for the above constructor.
     * A combination of default template parameters and auto-deduction guides works as expected,
     * independent of whether the second/optional parameter is specified or not, i.e. it is possible
     * to auto-deduct and overwrite a single template parameter out of the four if the optional parameter
     * is specified and use the default otherwise.
     */

    /*!\brief Construct from an existing stream and with specified format.
     * \tparam file_format   The format of the file in the stream, must satisfy seqan3::SequenceFileInputFormat.
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
    template <IStream2 stream_t,
              SequenceFileInputFormat file_format>
    sequence_file_input(stream_t                 & stream,
                        file_format        const & SEQAN3_DOXYGEN_ONLY(format_tag),
                        selected_field_ids const & SEQAN3_DOXYGEN_ONLY(fields_tag) = selected_field_ids{}) :
        primary_stream{&stream, stream_deleter_noop},
        format{detail::sequence_file_input_format<file_format>{}}
    {
        static_assert(meta::in<valid_formats, file_format>::value,
                      "You selected a format that is not in the valid_formats of this file.");

        // possibly add intermediate compression stream
        secondary_stream = detail::make_secondary_istream(*primary_stream);

        // buffer first record
        read_next_record();
    }

    //!\overload
    template <IStream2 stream_t,
              SequenceFileInputFormat file_format>
    sequence_file_input(stream_t                && stream,
                        file_format        const & SEQAN3_DOXYGEN_ONLY(format_tag),
                        selected_field_ids const & SEQAN3_DOXYGEN_ONLY(fields_tag) = selected_field_ids{}) :
        primary_stream{new stream_t{std::move(stream)}, stream_deleter_default},
        format{detail::sequence_file_input_format<file_format>{}}
    {
        static_assert(meta::in<valid_formats, file_format>::value,
                      "You selected a format that is not in the valid_formats of this file.");

        // possibly add intermediate compression stream
        secondary_stream = detail::make_secondary_istream(*primary_stream);

        // buffer first record
        read_next_record();
    }
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
     * \snippet test/snippet/io/sequence_file/sequence_file_input.cpp return_record
     *
     * It most situations using the iterator interface or a range-based for-loop are preferable to using front(),
     * because you can only move to the next record via the iterator.
     *
     * In any case, don't forget the reference! If you want to save the data from the record elsewhere, use move:
     *
     * \snippet test/snippet/io/sequence_file/sequence_file_input.cpp record_move
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
        return record_buffer;
    }
    //!\}

    /*!\name Tuple interface
     * \brief Provides functions for field-based ("column"-based) reading.
     * \{
     */
    //!\brief Read the entire file into internal buffers and retrieve the specified column.
    template <field f>
    friend auto & get(sequence_file_input & file)
    {
        static_assert(sequence_file_input::selected_field_ids::contains(f),
                      "You requested a field via get that was not selected for the file.");

        file.read_columns();

        return seqan3::get<f>(file.columns_buffer);
    }

    //!\copydoc get
    template <field f>
    friend auto && get(sequence_file_input && file)
    {
        return std::move(get<f>(file));
    }

    //!\copydoc get
    template <size_t i>
    friend auto & get(sequence_file_input & file)
    {
        static_assert(i < sequence_file_input::selected_field_ids::as_array.size(),
                      "You requested a field number larger than the number of selected fields for the file.");
        file.read_columns();

        return std::get<i>(file.columns_buffer);
    }

    //!\copydoc get
    template <size_t i>
    friend auto && get(sequence_file_input && file)
    {
        return std::move(get<i>(file));
    }

    //!\copydoc get
    template <typename t>
    friend auto & get(sequence_file_input & file)
    {
        file.read_columns();

        return std::get<t>(file.columns_buffer);
    }

    //!\copydoc get
    template <typename t>
    friend auto && get(sequence_file_input && file)
    {
        return std::move(get<t>(file));
    }
    //!\}

    //!\brief The options are public and its members can be set directly.
    sequence_file_input_options<typename traits_type::sequence_legal_alphabet,
                             selected_field_ids::contains(field::SEQ_QUAL)> options;

protected:
    //!\privatesection
    /*!\name Data buffers
     * \{
     */
    //!\brief Buffer for a single record.
    record_type record_buffer;
    //!\brief Buffer of the entire file in columns.
    file_as_tuple_type columns_buffer;
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

    //!\brief File is at position 1 behind the last record.
    bool at_end{false};

    //!\brief Type of the format, an std::variant over the `valid_formats`.
    using format_type = typename detail::variant_from_tags<valid_formats, detail::sequence_file_input_format>::type;
    //!\brief The actual std::variant holding a pointer to the detected/selected format.
    format_type format;
    //!\}

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

        assert(!format.valueless_by_exception());
        std::visit([&] (auto & f)
        {
            // read new record
            if constexpr (selected_field_ids::contains(field::SEQ_QUAL))
            {
                f.read(*secondary_stream,
                       options,
                       detail::get_or_ignore<field::SEQ_QUAL>(record_buffer),
                       detail::get_or_ignore<field::ID>(record_buffer),
                       detail::get_or_ignore<field::SEQ_QUAL>(record_buffer));
            }
            else
            {
                f.read(*secondary_stream,
                       options,
                       detail::get_or_ignore<field::SEQ>(record_buffer),
                       detail::get_or_ignore<field::ID>(record_buffer),
                       detail::get_or_ignore<field::QUAL>(record_buffer));
            }
        }, format);
    }

    //!\brief Read the entire file into the internal column buffers.
    void read_columns()
    {
        //TODO don't do multiple visits
        //TODO create specialised version for concatenated_sequences where we append on the concat
        auto & sequence_column_buffer = detail::get_or_ignore<field::SEQ>(columns_buffer);
        auto &       id_column_buffer = detail::get_or_ignore<field::ID>(columns_buffer);
        auto &     qual_column_buffer = detail::get_or_ignore<field::QUAL>(columns_buffer);
        auto & seq_qual_column_buffer = detail::get_or_ignore<field::SEQ_QUAL>(columns_buffer);

        // read the remaining records and split into column buffers
        for (auto & rec : *this)
        {
            if constexpr (selected_field_ids::contains(field::SEQ))
                sequence_column_buffer.push_back(std::move(seqan3::get<field::SEQ>(rec)));
            if constexpr (selected_field_ids::contains(field::ID))
                id_column_buffer.push_back(std::move(seqan3::get<field::ID>(rec)));
            if constexpr (selected_field_ids::contains(field::QUAL))
                qual_column_buffer.push_back(std::move(seqan3::get<field::QUAL>(rec)));
            if constexpr (selected_field_ids::contains(field::SEQ_QUAL))
                seq_qual_column_buffer.push_back(std::move(seqan3::get<field::SEQ_QUAL>(rec)));
        }
    }

    //!\brief Befriend iterator so it can access the buffers.
    friend iterator;
};

/*!\name Type deduction guides
 * \relates seqan3::sequence_file_input
 * \{
 */

//!\brief Deduction of the selected fields, the file format and the stream type.
template <IStream2                           stream_type,
          SequenceFileInputFormat            file_format,
          detail::Fields                     selected_field_ids>
sequence_file_input(stream_type && stream,
                    file_format const &,
                    selected_field_ids const &)
    -> sequence_file_input<typename sequence_file_input<>::traits_type,       // actually use the default
                           selected_field_ids,
                           type_list<file_format>,
                           typename std::remove_reference_t<stream_type>::char_type>;

//!\overload
template <IStream2                           stream_type,
          SequenceFileInputFormat            file_format,
          detail::Fields                     selected_field_ids>
sequence_file_input(stream_type & stream,
                    file_format const &,
                    selected_field_ids const &)
    -> sequence_file_input<typename sequence_file_input<>::traits_type,       // actually use the default
                           selected_field_ids,
                           type_list<file_format>,
                           typename std::remove_reference_t<stream_type>::char_type>;
//!\}

} // namespace seqan3

// ------------------------------------------------------------------
// std-overloads for the tuple-like interface
// ------------------------------------------------------------------

namespace std
{
/*!\brief Provides access to the number of elements in a tuple as a compile-time constant expression.
 * \implements seqan3::UnaryTypeTrait
 * \ingroup sequence
 * \see std::tuple_size_v
 */
template <seqan3::SequenceFileInputTraits                    traits_type,
          seqan3::detail::Fields                             selected_field_ids,
          seqan3::detail::TypeListOfSequenceFileInputFormats valid_formats,
          seqan3::Char                                       stream_char_t>
struct tuple_size<seqan3::sequence_file_input<traits_type, selected_field_ids, valid_formats, stream_char_t>>
{
    //!\brief The value equals the number of selected fields in the file.
    static constexpr size_t value = selected_field_ids::as_array.size();
};

/*!\brief Obtains the type of the specified element.
 * \implements seqan3::TransformationTrait
 * \ingroup sequence
 * \see [std::tuple_element](https://en.cppreference.com/w/cpp/utility/tuple/tuple_element)
 */
template <size_t                                             elem_no,
          seqan3::SequenceFileInputTraits                    traits_type,
          seqan3::detail::Fields                             selected_field_ids,
          seqan3::detail::TypeListOfSequenceFileInputFormats valid_formats,
          seqan3::Char                                       stream_char_t>
struct tuple_element<elem_no, seqan3::sequence_file_input<traits_type, selected_field_ids, valid_formats, stream_char_t>>
    : tuple_element<elem_no, typename seqan3::sequence_file_input<traits_type,
                                                               selected_field_ids,
                                                               valid_formats,
                                                               stream_char_t>::file_as_tuple_type>
{};

} // namespace std

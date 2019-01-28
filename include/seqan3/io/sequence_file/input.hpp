// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
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
#include <seqan3/core/metafunction/basic.hpp>
#include <seqan3/io/stream/concept.hpp>
#include <seqan3/io/exception.hpp>
#include <seqan3/io/filesystem.hpp>
#include <seqan3/io/record.hpp>
#include <seqan3/io/detail/in_file_iterator.hpp>
#include <seqan3/io/detail/misc_input.hpp>
#include <seqan3/io/detail/record.hpp>
#include <seqan3/io/sequence_file/input_format_concept.hpp>
#include <seqan3/io/sequence_file/format_fasta.hpp>
#include <seqan3/io/sequence_file/format_fastq.hpp>
#include <seqan3/range/container/concatenated_sequences.hpp>

namespace seqan3
{

// ----------------------------------------------------------------------------
// sequence_file_input_traits_concept
// ----------------------------------------------------------------------------

/*!\interface seqan3::sequence_file_input_traits_concept <>
 * \brief The requirements a traits_type for seqan3::sequence_file_input must meet.
 * \ingroup sequence
 */
/*!\name Requirements for seqan3::sequence_file_input_traits_concept
 * \brief You can expect these **member types** of all types that satisfy seqan3::sequence_file_input_traits_concept.
 * \memberof seqan3::sequence_file_input_traits_concept
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
 * \memberof seqan3::sequence_file_input_traits_concept
 * \brief Alphabet of the characters for the seqan3::field::SEQ; must satisfy seqan3::alphabet_concept.
 */
/*!\typedef using sequence_legal_alphabet
 * \memberof seqan3::sequence_file_input_traits_concept
 * \brief Intermediate alphabet for seqan3::field::SEQ; must satisfy seqan3::alphabet_concept and be convertible to
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
 * \memberof seqan3::sequence_file_input_traits_concept
 * \brief Type template of the seqan3::field::SEQ, a container template over `sequence_alphabet`;
 * must satisfy seqan3::sequence_container_concept.
 */
/*!\typedef using sequence_container_container
 * \memberof seqan3::sequence_file_input_traits_concept
 * \brief Type template of a column of seqan3::field::SEQ, a container template that can hold multiple
 * `sequence_container`; must satisfy seqan3::sequence_container_concept.
 */
/*!\typedef using id_alphabet
 * \memberof seqan3::sequence_file_input_traits_concept
 * \brief Alphabet of the characters for the seqan3::field::ID; must satisfy seqan3::alphabet_concept.
 */
/*!\typedef using id_container
 * \memberof seqan3::sequence_file_input_traits_concept
 * \brief Type template of the seqan3::field::ID, a container template over `id_alphabet`;
 * must satisfy seqan3::sequence_container_concept.
 */
/*!\typedef using id_container_container
 * \memberof seqan3::sequence_file_input_traits_concept
 * \brief Type template of a column of seqan3::field::ID, a container template that can hold multiple
 * `id_container`; must satisfy seqan3::sequence_container_concept.
 */
/*!\typedef using quality_alphabet
 * \memberof seqan3::sequence_file_input_traits_concept
 * \brief Alphabet of the characters for the seqan3::field::QUAL; must satisfy seqan3::quality_concept.
 */
/*!\typedef using quality_container
 * \memberof seqan3::sequence_file_input_traits_concept
 * \brief Type template of the seqan3::field::QUAL, a container template over `quality_alphabet`;
 * must satisfy seqan3::sequence_container_concept.
 */
/*!\typedef using quality_container_container
 * \memberof seqan3::sequence_file_input_traits_concept
 * \brief Type template of a column of seqan3::field::QUAL, a container template that can hold multiple
 * `quality_container`; must satisfy seqan3::sequence_container_concept.
 */
//!\}
//!\cond
template <typename t>
SEQAN3_CONCEPT sequence_file_input_traits_concept = requires (t v)
{
    requires alphabet_concept<typename t::sequence_alphabet>;
    requires alphabet_concept<typename t::sequence_legal_alphabet>;
    requires explicitly_convertible_to_concept<typename t::sequence_legal_alphabet, typename t::sequence_alphabet>;
    requires sequence_container_concept<typename t::template sequence_container<typename t::sequence_alphabet>>;
    requires sequence_container_concept<typename t::template sequence_container_container<
        typename t::template sequence_container<typename t::sequence_alphabet>>>;

    requires alphabet_concept<typename t::id_alphabet>;
    requires sequence_container_concept<typename t::template id_container<typename t::id_alphabet>>;
    requires sequence_container_concept<typename t::template id_container_container<typename t::template id_container<
        typename t::id_alphabet>>>;

    requires quality_concept<typename t::quality_alphabet>;
    requires sequence_container_concept<typename t::template quality_container<typename t::quality_alphabet>>;
    requires sequence_container_concept<typename t::template quality_container_container<
        typename t::template quality_container<typename t::quality_alphabet>>>;
};
//!\endcond

// ----------------------------------------------------------------------------
// sequence_file_input_default_traits
// ----------------------------------------------------------------------------

/*!\brief The default traits for seqan3::sequence_file_input
 * \implements sequence_file_input_traits_concept
 * \ingroup sequence
 *
 * \details
 *
 * If you wish to change a single or a few types from the default, just inherit from this class and
 * "overwrite" the respective type definitions.
 *
 * This example will make the file read into a smaller alphabet and a compressed container:
 *
 * ```cpp
 * struct my_traits : sequence_file_input_default_traits_dna
 * {
 *     using sequence_alphabet = dna4;                        // instead of dna5
 *
 *     template <typename alph>
 *     using sequence_container = bitcompressed_vector<alph>; // must be defined as a template!
 * };
 *
 * sequence_file_input<my_traits> fin{"/tmp/my.fasta"};
 *
 * //...
 * ```
 */
struct sequence_file_input_default_traits_dna
{
    /*!\name Member types
     * \brief Definitions to satisfy seqan3::sequence_file_input_traits_concept.
     * \{
     */
    using sequence_alphabet                 = dna5;
    using sequence_legal_alphabet           = dna15;
    template <typename _sequence_alphabet>
    using sequence_container                = std::vector<_sequence_alphabet>;
    template <typename _sequence_container>
    using sequence_container_container      = concatenated_sequences<_sequence_container>;

    using id_alphabet                       = char;
    template <typename _id_alphabet>
    using id_container                      = std::basic_string<_id_alphabet>;
    template <typename _id_container>
    using id_container_container            = concatenated_sequences<_id_container>;

    using quality_alphabet                  = phred42;
    template <typename _quality_alphabet>
    using quality_container                 = std::vector<_quality_alphabet>;
    template <typename _quality_container>
    using quality_container_container       = concatenated_sequences<_quality_container>;
    //!\}
};

//!\brief A traits type that specifies input as amino acids.
//!\ingroup sequence
struct sequence_file_input_default_traits_aa : sequence_file_input_default_traits_dna
{
    /*!\name Member types
     * \brief Definitions to satisfy seqan3::sequence_file_input_traits_concept.
     * \{
     */
    using sequence_alphabet = aa27;
    using sequence_legal_alphabet = aa27;
    //!\}
};

// ----------------------------------------------------------------------------
// sequence_file_input
// ----------------------------------------------------------------------------

/*!\brief A class for reading sequence files, e.g. FASTA, FASTQ ...
 * \ingroup sequence
 * \tparam traits_type          An auxiliary type that defines certain member types and constants, must satisfy
 * seqan3::sequence_file_input_traits_concept.
 * \tparam selected_field_ids   A seqan3::fields type with the list and order of desired record entries; all fields
 * must be in seqan3::sequence_file_input::field_ids.
 * \tparam valid_formats        A seqan3::type_list of the selectable formats (each must meet
 * seqan3::sequence_file_input_format_concept).
 * \tparam stream_char_type     The type of the underlying stream device(s); must model seqan3::char_concept.
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
 * ```cpp
 * sequence_file_input fin{"/tmp/my.fasta"}; // FastA with DNA sequences assumed, regular std::ifstream taken as stream
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
 * sequence_file_input fin{std::move(iss), sequence_file_format_fasta{}};
 * //              ^ no need to specify the template arguments
 * ```
 *
 * Note that this is not the same as writing `sequence_file_input<>` (with angle brackets). In the latter case they are
 * explicitly set to their default values, in the former case
 * [automatic deduction](http://en.cppreference.com/w/cpp/language/class_template_argument_deduction) happens which
 * chooses different parameters depending on the constructor arguments. For opening from file, `sequence_file_input<>`
 * would have also worked, but for opening from stream it would not have.
 *
 * In some cases, you do need to specify the arguments, e.g. if you want to read amino acids:
 *
 * ```cpp
 * sequence_file_input<sequence_file_default_traits_aa> fin{"/tmp/my.fasta"};
 * ```
 *
 * You can define your own traits type to further customise the types used by and returned by this class, see
 * seqan3::sequence_file_default_traits_dna for more details. As mentioned above, specifying at least one
 * template parameter yourself means that you loose automatic deduction so if you want to read amino acids **and**
 * want to read from a string stream you need to give all types yourself:
 *
 * ```cpp
 *
 *  // ... input had amino acid sequences
 * std::istringstream iss(input);
 *
 * sequence_file_input<sequence_file_default_traits_aa,
 *                  fields<field::SEQ, field::ID, field::QUAL>,
 *                  type_list<sequence_file_format_fasta>,
 *                  std::istringstream> fin{std::move(iss), sequence_file_format_fasta{}};
 * ```
 *
 * ### Reading record-wise
 *
 * You can iterate over this file record-wise:
 *
 * ```cpp
 * sequence_file_input fin{"/tmp/my.fasta"};
 *
 * for (auto & rec : fin)
 * {
 *     std::cout << "ID:  " << get<field::ID>(rec) << '\n';
 *     std::cout << "SEQ: " << (get<field::SEQ>(rec) | view::to_char) << '\n'; // sequence is converted to char on-the-fly
 *     // a quality field also exists, but is not printed, because we know it's empty for FastA files.
 * }
 * ```
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
 * ```cpp
 * sequence_file_input fin{"/tmp/my.fasta"};
 *
 * using record_type = typename decltype(fin)::record_type;
 * std::vector<record_type> records;
 *
 * for (auto & rec : fin)
 *     records.push_back(std::move(rec));
 * ```
 *
 * ### Reading record-wise (decomposed records)
 *
 * Instead of using `get` on the record, you can also use
 * [structured bindings](http://en.cppreference.com/w/cpp/language/structured_binding)
 * to decompose the record into its elements:
 *
 * ```cpp
 * sequence_file_input fin{"/tmp/my.fasta"};
 *
 * for (auto & [ seq, id, qual ] : fin)
 * {
 *     std::cout << "ID:  " << id << '\n';
 *     std::cout << "SEQ: " << (seq | view::to_char) << '\n'; // sequence is converted to char on-the-fly
 *     // qual is empty for FastA files
 * }
 * ```
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
 * ```cpp
 * sequence_file_input fin{"/tmp/my.fasta", fields<field::ID, field::SEQ_QUAL>{}};
 *
 * for (auto & [ id, seq_qual ] : fin) // note that the order is now different, "id" comes first, because it was specified first
 * {
 *     std::cout << "ID:  " << id << '\n';
 *     // sequence and qualities are part of the same vector, of type std::vector<dna5q>
 *     std::cout << "SEQ: "  << (seq | view::get<0> | view::to_char) << '\n'; // sequence string is extracted and converted to char  on-the-fly
 *     std::cout << "QUAL: " << (seq | view::get<1> | view::to_char) << '\n'; // quality string is extracted and converted to char  on-the-fly
 * }
 * ```
 *
 * When reading a file, all fields not present in the file (but requested implicitly or via the `selected_field_ids`
 * parameter) are ignored.
 *
 * ### Views on files
 *
 * Since SeqAn files are ranges, you can also create views over files. A useful example is to filter the records
 * based on certain criteria, e.g. minimum length of the sequence field:
 *
 * ```cpp
 * sequence_file_input fin{"/tmp/my.fasta"};
 *
 * auto minimum_length5_filter = view::filter([] (auto const & rec)
 * {
 *     return size(get<field::SEQ>(rec)) >= 5;
 * });
 *
 * for (auto & rec : fin | minimum_length5_filter) // only record with sequence length >= 5 will "appear"
 * {
 *     // ...
 * }
 * ```
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
 * ```cpp
 *
 * struct data_storage_t
 * {
 *     concatenated_sequences<dna5_vector>  sequences;
 *     concatenated_sequences<std::string>  ids;
 * };
 *
 * data_storage_t data_storage; // a global or globally used variable in your program
 *
 * // ... in your file reading function:
 *
 * sequence_file_input fin{"/tmp/my.fasta"};
 *
 * data_storage.sequences = std::move(get<field::SEQ>(fin)); // we move the buffer directly into our storage
 * data_storage.ids = std::move(get<field::ID>(fin)); // we move the buffer directly into our storage
 * ```
 *
 * Note that for this to make sense, your storage data types need to be identical to the corresponding column types
 * of the file. If you require different column types you can specify you own traits, see
 * seqan3::sequence_file_input_traits_concept.
 *
 * ### Formats
 *
 * TODO give overview of formats, once they are all implemented
 */

template <
    sequence_file_input_traits_concept                       traits_type_        = sequence_file_input_default_traits_dna,
    detail::fields_concept                                   selected_field_ids_ = fields<field::SEQ,
                                                                                          field::ID,
                                                                                          field::QUAL>,
    detail::type_list_of_sequence_file_input_formats_concept valid_formats_      = type_list<sequence_file_format_fasta,
                                                                                          sequence_file_format_fastq
                                                                                         /*, ...*/>,
    char_concept                                             stream_char_type_   = char>
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
    using sequence_quality_column_type  = typename traits_type::template id_container_container<sequence_quality_type>;
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
    using sentinel          = std::ranges::default_sentinel;
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
    sequence_file_input(filesystem::path filename,
                        selected_field_ids const & SEQAN3_DOXYGEN_ONLY(fields_tag) = selected_field_ids{}) :
        primary_stream{new std::ifstream{filename, std::ios_base::in | std::ios::binary}, stream_deleter_default}
    {
        if (!primary_stream->good())
            throw file_open_error{"Could not open file for reading."};

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
     * \tparam file_format   The format of the file in the stream, must satisfy seqan3::sequence_file_input_format_concept.
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
    template <istream_concept2 stream_t,
              sequence_file_input_format_concept file_format>
    sequence_file_input(stream_t                 & stream,
                        file_format        const & SEQAN3_DOXYGEN_ONLY(format_tag),
                        selected_field_ids const & SEQAN3_DOXYGEN_ONLY(fields_tag) = selected_field_ids{}) :
        primary_stream{&stream, stream_deleter_noop},
        format{file_format{}}
    {
        static_assert(meta::in<valid_formats, file_format>::value,
                      "You selected a format that is not in the valid_formats of this file.");

        // possibly add intermediate compression stream
        secondary_stream = detail::make_secondary_istream(*primary_stream);

        // buffer first record
        read_next_record();
    }

    //!\overload
    template <istream_concept2 stream_t,
              sequence_file_input_format_concept file_format>
    sequence_file_input(stream_t                && stream,
                        file_format        const & SEQAN3_DOXYGEN_ONLY(format_tag),
                        selected_field_ids const & SEQAN3_DOXYGEN_ONLY(fields_tag) = selected_field_ids{}) :
        primary_stream{new stream_t{std::move(stream)}, stream_deleter_default},
        format{file_format{}}
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
     * ```cpp
     * sequence_file_input fin{"/tmp/my.fasta"};
     * auto it = begin(fin);
     *
     * // the following are equivalent:
     * auto & rec0 = *it;
     * auto & rec1 = fin.front();
     *
     * // both become invalid after incrementing "it"!
     * ```
     *
     * It most situations using the iterator interface or a range-based for-loop are preferable to using front(),
     * because you can only move to the next record via the iterator.
     *
     * In any case, don't forget the reference! If you want to save the data from the record elsewhere, use move:
     *
     * ```cpp
     * sequence_file_input fin{"/tmp/my.fasta"};
     *
     * auto rec0 = std::move(fin.front());
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
    using format_type = detail::transfer_template_args_onto_t<valid_formats, std::variant>;
    //!\brief The actual std::variant holding a pointer to the detected/selected format.
    format_type format;
    //!\}

    //!\brief Tell the format to move to the next record and update the buffer.
    void read_next_record()
    {
        if (at_end)
            return;

        // clear the record
        record_buffer.clear();

        // at end if we could not read further
        if (secondary_stream->eof())
        {
            at_end = true;
            return;
        }

        assert(!format.valueless_by_exception());
        std::visit([&] (sequence_file_input_format_concept & f)
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
template <istream_concept2                   stream_type,
          sequence_file_input_format_concept file_format,
          detail::fields_concept             selected_field_ids>
sequence_file_input(stream_type && stream,
                    file_format const &,
                    selected_field_ids const &)
    -> sequence_file_input<typename sequence_file_input<>::traits_type,       // actually use the default
                           selected_field_ids,
                           type_list<file_format>,
                           typename std::remove_reference_t<stream_type>::char_type>;

template <istream_concept2                   stream_type,
          sequence_file_input_format_concept file_format,
          detail::fields_concept             selected_field_ids>
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
//!\brief std::tuple_size overload for column-like access. [metafunction specialisation for seqan3::sequence_file_input]
template <seqan3::sequence_file_input_traits_concept                       traits_type,
          seqan3::detail::fields_concept                                selected_field_ids,
          seqan3::detail::type_list_of_sequence_file_input_formats_concept valid_formats,
          seqan3::char_concept                                             stream_char_t>
struct tuple_size<seqan3::sequence_file_input<traits_type, selected_field_ids, valid_formats, stream_char_t>>
{
    //!\brief The value equals the number of selected fields in the file.
    static constexpr size_t value = selected_field_ids::as_array.size();
};

//!\brief std::tuple_element overload for column-like access. [metafunction specialisation for seqan3::sequence_file_input]
template <size_t                                                        elem_no,
          seqan3::sequence_file_input_traits_concept                       traits_type,
          seqan3::detail::fields_concept                                selected_field_ids,
          seqan3::detail::type_list_of_sequence_file_input_formats_concept valid_formats,
          seqan3::char_concept                                             stream_char_t>
struct tuple_element<elem_no, seqan3::sequence_file_input<traits_type, selected_field_ids, valid_formats, stream_char_t>>
    : tuple_element<elem_no, typename seqan3::sequence_file_input<traits_type,
                                                               selected_field_ids,
                                                               valid_formats,
                                                               stream_char_t>::file_as_tuple_type>
{};

} // namespace std

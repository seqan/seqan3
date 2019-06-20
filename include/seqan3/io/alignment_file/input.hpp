// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::alignment_file_input and corresponding traits classes.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <cassert>
#include <fstream>
#include <string>
#include <variant>
#include <vector>

#include <seqan3/alphabet/adaptation/char.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/alphabet/quality/qualified.hpp>
#include <seqan3/core/concept/tuple.hpp>
#include <seqan3/core/type_traits/basic.hpp>
#include <seqan3/core/type_traits/transformation_trait_or.hpp>
#include <seqan3/io/alignment_file/input_format_concept.hpp>
#include <seqan3/io/alignment_file/format_bam.hpp>
#include <seqan3/io/alignment_file/format_sam.hpp>
#include <seqan3/io/alignment_file/misc.hpp>
#include <seqan3/io/detail/in_file_iterator.hpp>
#include <seqan3/io/detail/misc_input.hpp>
#include <seqan3/io/detail/record.hpp>
#include <seqan3/io/exception.hpp>
#include <seqan3/io/record.hpp>
#include <seqan3/io/stream/concept.hpp>
#include <seqan3/range/container/concatenated_sequences.hpp>
#include <seqan3/range/decorator/gap_decorator.hpp>
#include <seqan3/range/view/repeat_n.hpp>
#include <seqan3/range/view/slice.hpp>
#include <seqan3/std/concepts>
#include <seqan3/std/filesystem>
#include <seqan3/std/ranges>

namespace seqan3
{

// ---------------------------------------------------------------------------------------------------------------------
// AlignmentFileInputTraits
// ---------------------------------------------------------------------------------------------------------------------

/*!\interface seqan3::AlignmentFileInputTraits <>
 * \brief The requirements a traits_type for seqan3::alignment_file_input must meet.
 * \ingroup alignment_file
 */
/*!\name Requirements for seqan3::AlignmentFileInputTraits
 * \brief You can expect these **member types** of all types that model seqan3::AlignmentFileInputTraits.
 * \memberof seqan3::AlignmentFileInputTraits
 * \{
 */
/*!\typedef using sequence_alphabet
 * \brief Alphabet of the characters for the seqan3::field::SEQ; must model seqan3::Alphabet.
 */
/*!\typedef using sequence_legal_alphabet
 * \brief Intermediate alphabet for seqan3::field::SEQ; must model seqan3::Alphabet and be convertible to
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
 * must model seqan3::SequenceContainer.
 */
/*!\typedef using id_container
 * \brief Type template of the seqan3::field::ID, a container template over `char`;
 * must model seqan3::SequenceContainer.
 */
/*!\typedef using quality_alphabet
 * \brief Alphabet of the characters for the seqan3::field::QUAL; must model seqan3::WritableQualityAlphabet.
 */
/*!\typedef using quality_container
 * \brief Type template of the seqan3::field::QUAL, a container template over `quality_alphabet`;
 * must model seqan3::SequenceContainer.
 */
/*!\typedef using ref_sequences
 * \brief The type of range over reference sequences; must model std::ranges::ForwardRange,
 *        the value_type must also model std::ranges::ForwardRange, and the value type of the value type
 *        must model seqan3::Alphabet (e.g. std::vector<std::vector<dna4>>).
 *
 * \attention This type is the first template parameter of the alignment_file_default_traits and should not be manually
 *            configured in order to allow for automatic type deduction from reference information input on
 *            construction.
 */
/*!\typedef using ref_ids
 * \brief The type of range over reference sequences; must model std::ranges::ForwardRange,
 *        the value_type must also model std::ranges::ForwardRange, and the value type of the value type
 *        must model seqan3::Alphabet (e.g. std::vector<string>).
 *
 * \attention This type is the second template parameter of the alignment_file_default_traits and should not be manually
 *            configured in order to allow for automatic type deduction from reference information input on
 *            construction.
 */
//!\}
//!\cond
template <typename t>
SEQAN3_CONCEPT AlignmentFileInputTraits = requires (t v)
{
    // field::SEQ
    requires WritableAlphabet<typename t::sequence_alphabet>;
    requires WritableAlphabet<typename t::sequence_legal_alphabet>;
    requires ExplicitlyConvertibleTo<typename t::sequence_legal_alphabet, typename t::sequence_alphabet>;
    requires SequenceContainer<typename t::template sequence_container<typename t::sequence_alphabet>>;

    // field::ID
    requires SequenceContainer<typename t::template id_container<char>>;

    // field::QUAL
    requires WritableQualityAlphabet<typename t::quality_alphabet>;
    requires SequenceContainer<typename t::template quality_container<typename t::quality_alphabet>>;

    // field::REF_SEQ
    // either ref_info_not_given or a range over ranges over Alphabet (e.g. std::vector<dna4_vector>)
    requires std::Same<typename t::ref_sequences, ref_info_not_given> ||
         (std::ranges::ForwardRange<typename t::ref_sequences> &&
         std::ranges::ForwardRange<detail::transformation_trait_or_t<reference<typename t::ref_sequences>, dna4_vector>> &&
         Alphabet<reference_t<detail::transformation_trait_or_t<reference<typename t::ref_sequences>, dna4_vector>>>);

    // field::REF_ID
    requires Alphabet<reference_t<reference_t<typename t::ref_ids>>> &&
             (!std::Same<typename t::ref_sequences, ref_info_not_given> ||
              WritableAlphabet<reference_t<reference_t<typename t::ref_ids>>>);
    requires std::ranges::ForwardRange<reference_t<typename t::ref_ids>>;
    requires std::ranges::ForwardRange<typename t::ref_ids>;

    // field::OFFSET is fixed to int32_t
    // field::REF_OFFSET is fixed to std::optional<int32_t>
    // field::FLAG is fixed to uint16_t
    // field::MAPQ is fixed to uint8_t
    // field::EVALUE is fixed to double
    // field::BITSCORE is fixed to double
    // field::MATE is fixed to std::tuple<ref_id_container<ref_id_alphabet>, ref_offset_type, int32_t>

    // field::ALIGNMENT
    // the alignment type cannot be configured.
    // Type of tuple entry 1 (reference) is set to
    // 1) a std::ranges::subrange over value_type_t<typename t::ref_sequences> if reference information was given
    // or 2) a "dummy" sequence type:
    // view::repeat_n(sequence_alphabet{}, size_t{}) | std::view::transform(detail::access_restrictor_fn{})
    // Type of tuple entry 2 (query) is set to
    // 1) a std::ranges::subrange over value_type_t<typename t::ref_sequences> if reference information was given
    // or 2) a "dummy" sequence type:
};
//!\endcond

// ---------------------------------------------------------------------------------------------------------------------
// alignment_file_input_default_traits
// ---------------------------------------------------------------------------------------------------------------------

/*!\brief The default traits for seqan3::alignment_file_input
 * \implements AlignmentFileInputTraits
 * \ingroup alignment_file
 * \tparam ref_sequences_t A range over reference sequences. This type is automatically deduced on construction.
 * \tparam ref_ids_t       A range over reference ids. This type is automatically deduced on construction.
 *
 * \details
 *
 * If you wish to change a single or a few types from the default, just inherit from this class and
 * "overwrite" the respective type definitions.
 *
 * This example will make the file read into a smaller alphabet and a compressed container:
 *
 * \snippet test/snippet/io/alignment_file/alignment_file_input.cpp my_traits
 */
template <typename ref_sequences_t = ref_info_not_given, typename ref_ids_t = std::deque<std::string>>
struct alignment_file_input_default_traits
{
    /*!\name Member types
     * \brief Definitions to model seqan3::AlignmentFileInputTraits.
     * \{
     */

    //!\brief The sequence alphabet is seqan3::dna5.
    using sequence_alphabet                     = dna5;

    //!\brief The legal sequence alphabet for parsing is seqan3::dna15.
    using sequence_legal_alphabet               = dna15;

    //!\brief The container for a sequence is std::vector.
    template <typename _sequence_alphabet>
    using sequence_container                    = std::vector<_sequence_alphabet>;

    //!\brief The string type for an identifier is std::basic_string.
    template <typename _id_alphabet>
    using id_container                          = std::basic_string<_id_alphabet>;

    //!\brief The alphabet for a quality annotation is seqan3::phred42.
    using quality_alphabet                      = phred42;

    //!\brief The string type for a quality annotation is std::vector.
    template <typename _quality_alphabet>
    using quality_container                     = std::vector<_quality_alphabet>;

    //!\brief The type of the reference sequences is deduced on construction.
    using ref_sequences                         = ref_sequences_t;

    //!\brief The type of the reference identifiers is deduced on construction.
    using ref_ids                               = ref_ids_t;
    //!\}
};

// ---------------------------------------------------------------------------------------------------------------------
// alignment_file_input
// ---------------------------------------------------------------------------------------------------------------------

/*!\brief A class for reading alignment files, e.g. SAM, BAM, BLAST ...
 * \ingroup alignment_file
 * \tparam traits_type          An auxiliary type that defines certain member types and constants, must model
 *                              seqan3::AlignmentFileInputTraits.
 * \tparam selected_field_ids   A seqan3::fields type with the list and order of desired record entries; all fields
 *                              must be in seqan3::alignment_file_input::field_ids.
 * \tparam valid_formats        A seqan3::type_list of the selectable formats (each must meet
 *                              seqan3::AlignmentFileInputFormat).
 * \tparam stream_char_type     The type of the underlying stream device(s); must model std::Integral.
 *
 * \details
 *
 * ### Introduction
 *
 * Alignment files are primarily used to store pairwise alignments of two biological sequences and often come with
 * many additional information. Well-known formats include the SAM/BAM format used to store read mapping data or the
 * BLAST format that stores the results of a query search against a data base.
 *
 * The Alignment file abstraction supports reading 14 different fields:
 *
 *   1. seqan3::field::SEQ
 *   2. seqan3::field::ID
 *   3. seqan3::field::OFFSET
 *   4. seqan3::field::REF_SEQ
 *   5. seqan3::field::REF_ID
 *   6. seqan3::field::REF_OFFSET
 *   7. seqan3::field::ALIGNMENT
 *   8. seqan3::field::MAPQ
 *   9. seqan3::field::QUAL
 *   10. seqan3::field::FLAG
 *   11. seqan3::field::MATE
 *   12. seqan3::field::TAGS
 *   13. seqan3::field::EVALUE
 *   14. seqan3::field::BIT_SCORE
 *
 * There exists one more field for alignment files, the seqan3::field::HEADER_PTR, but this field is mostly used
 * internally. Please see the seqan3::alignment_file_output::header member function for details on how to access the
 * seqan3::alignment_file_header of the file.)
 *
 * All of these fields are retrieved by default (and in that order).
 * Note that some of the fields are specific to the SAM format (e.g. seqan3::field::FLAG) while others are specific to
 * BLAST format (e.g. seqan3::field::BIT_SCORE). Please see the corresponding formats for more details
 * (seqan3::format_sam).
 *
 * ### Construction and specialisation
 *
 * This class comes with four constructors: One for construction from a file name, one for construction from
 * an existing stream and a known format and both of the former with or without additional reference information.
 *
 * Constructing from a file name automatically picks the format based on the extension
 * of the file name. Constructing from a stream can be used if you have a non-file stream, like std::cin or
 * std::istringstream, that you want to read from and/or if you cannot use file-extension based detection,
 * but know that your input file has a certain format.
 *
 * The reference information is specific to the SAM format. The SAM format only stores a "semi-alignment" meaning that
 * it has the query sequence and the cigar string representing the gap information but not the reference information.
 * If you want to retrieve valid/full alignments, you need to pass the corresponding reference information:
 *
 * - ref_ids: The name of the references, e.g. "chr1", "chr2", ...
 * - ref_sequences: The reference sequence information **in the same order as the ref_ids**.
 *
 * In most cases the template parameters are deduced automatically:
 *
 * \snippet test/snippet/io/alignment_file/alignment_file_input.cpp construction_from_filename
 *
 * Reading from an std::istringstream:
 *
 * \snippet test/snippet/io/alignment_file/alignment_file_input.cpp construction_from_stream
 *
 * Note that this is not the same as writing `alignment_file_input<>` (with angle brackets). In the latter case they
 * are explicitly set to their default values, in the former case
 * [automatic deduction](http://en.cppreference.com/w/cpp/language/class_template_argument_deduction) happens which
 * chooses different parameters depending on the constructor arguments. For opening from file, `alignment_file_input<>`
 * would have also worked, but for opening from stream it would not have.
 *
 * You can define your own traits type to further customise the types used by and returned by this class, see
 * seqan3::alignment_file_input_default_traits for more details. As mentioned above, specifying at least one
 * template parameter yourself means that you loose automatic deduction. The following is equivalent to the automatic
 * type deduction example with a stream from above:
 *
 * \snippet test/snippet/io/alignment_file/alignment_file_input.cpp construction_without_automatic_type_deduction
 *
 * ### Reading record-wise
 *
 * You can iterate over this file record-wise:
 *
 * \snippet test/snippet/io/alignment_file/alignment_file_input.cpp reading_range_based_for_loop
 *
 * In the above example, `rec` has the type \ref record_type which is a specialisation of seqan3::record and behaves
 * like an std::tuple (that's why we can access it via `get`). Instead of using the seqan3::field based interface on
 * the record, you could also use `std::get<0>` or even `std::get<dna4_vector>` to retrieve the sequence, but it is
 * not recommended, because it is more error-prone.
 *
 * *Note:* It is important to write `auto &` and not just `auto`, otherwise you will copy the record on every iteration.
 * Since the buffer gets "refilled" on every iteration, you can also move the data out of the record if you want
 * to store it somewhere without copying:
 *
 * \snippet test/snippet/io/alignment_file/alignment_file_input.cpp reading_move_record
 *
 * ### Reading record-wise (custom fields)
 *
 * If you want to skip specific fields from the record you can pass a non-empty fields trait object to the
 * seqan3::alignment_file_input constructor to select the fields that should be read from the input. For example,
 * you may only be interested in the mapping flag and mapping quality of your SAM data to get some statistics.
 * The following snippets demonstrate the usage of such a fields trait object.
 *
 * \snippet test/snippet/io/alignment_file/alignment_file_input.cpp reading_custom_fields
 *
 * When reading a file, all fields not present in the file (but requested implicitly or via the `selected_field_ids`
 * parameter) are ignored and the respective value in the record stays empty.
 *
 * ### Reading record-wise (decomposed records)
 *
 * Instead of using `get` on the record, you can also use
 * [structured bindings](http://en.cppreference.com/w/cpp/language/structured_binding)
 * to decompose the record into its elements. Considering the example of reading only the flag and mapping quality
 * like before you can also write:
 *
 * \snippet test/snippet/io/alignment_file/alignment_file_input.cpp reading_structured_bindings
 *
 * In this case you immediately get the two elements of the tuple: `flag` of \ref flag_type and `mapq` of
 * \ref mapq_type. **But beware: with structured bindings you do need to get the order of elements correctly!**
 *
 * ### Views on files
 *
 * Since SeqAn files are ranges, you can also create views over files. A useful example is to filter the records
 * based on certain criteria, e.g. minimum length of the sequence field:
 *
 * \snippet test/snippet/io/alignment_file/alignment_file_input.cpp reading_filter
 *
 * ### End of file
 *
 * You can check whether a file is at its end by comparing begin() and end() (if they are the same, the file is
 * at its end).
 *
 * ### Formats
 *
 * We currently support reading the following formats:
 *   * seqan3::format_sam
 *   * seqan3::format_bam
 */
template <
    AlignmentFileInputTraits                     traits_type_        = alignment_file_input_default_traits<>,
    detail::Fields                               selected_field_ids_ = fields<field::SEQ,
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
    detail::TypeListOfAlignmentFileInputFormats  valid_formats_    = type_list<format_sam, format_bam>,
    std::Integral                                stream_char_type_ = char>
class alignment_file_input
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

private:
    //!\brief The dummy ref sequence type if no reference information were given.
    using dummy_ref_type = decltype(view::repeat_n(typename traits_type::sequence_alphabet{}, size_t{}) |
                                    std::view::transform(detail::access_restrictor_fn{}));
public:
    /*!\name Field types and record type
     * \brief These types are relevant for record/row-based reading; they may be manipulated via the \ref traits_type
     * to achieve different storage behaviour.
     * \{
     */
    //!\brief The type of field::SEQ (default std::vector<seqan3::dna5>).
    using sequence_type            = typename traits_type::template sequence_container<
                                         typename traits_type::sequence_alphabet>;
    //!\brief The type of field::ID (default std::string by default).
    using id_type                  = typename traits_type::template id_container<char>;
    //!\brief The type of field::OFFSET is fixed to int32_t.
    using offset_type              = int32_t;
    /*!\brief The type of field::REF_SEQ (default depends on construction).
     *
     * If no reference information are given on construction, this type deduces to a sized view that throws on
     * access (since there is nothing to access anyway). If the reference information are given, the type is deduced
     * to a view over the given input reference sequence type such that no sequence information is copied.
     */
    using ref_sequence_type = std::conditional_t<std::Same<typename traits_type::ref_sequences, ref_info_not_given>,
                                                 dummy_ref_type,
                                                 decltype(std::declval<
                                                     detail::transformation_trait_or_t<
                                                        seqan3::reference<typename traits_type::ref_sequences const>,
                                                        dummy_ref_type> /* does not matter as type is not chosen */
                                                     >() | view::slice(0, 0))>;
    /*!\brief The type of field::REF_ID is fixed to std::optional<int32_t>.
     *
     * To be consistent with the BAM format, the field::REF_ID will hold the index to the actual reference
     * information stored in the header. If a read is unmapped, the optional will remain valueless.
     *
     * \attention SeqaAn3 transforms the 1-based SAM format position into a 0-based position.
     */
    using ref_id_type              = std::optional<int32_t>;
    /*!\brief The type of field::REF_OFFSET is fixed to an std::optional<int32_t>.
     *
     * The SAM format is 1-based and a 0 in the ref_offset field indicated an unmapped read. Since we convert 1-based
     * positions to 0-based positions when reading the SAM format, we model the ref_offset_type as an std::optional.
     * If the input value is 0, the std::optional will remain valueless.
     */
    using ref_offset_type          = std::optional<int32_t>;
    //!\brief The type of field::MAPQ is fixed to uint8_t.
    using mapq_type                = uint8_t;
    //!\brief The type of field::QUAL (default std::vector<seqan3::phred42>).
    using quality_type             = typename traits_type::template quality_container<
                                         typename traits_type::quality_alphabet>;
    //!\brief The type of field::FLAG is fixed to uint16_t.
    using flag_type                = uint16_t;
    //!\brief The type of field::MATE is fixed to std::tuple<ref_id_type, ref_offset_type, int32_t>).
    using mate_type                = std::tuple<ref_id_type, ref_offset_type, int32_t>;
    //!\brief The type of field::EVALUE is fixed to double.
    using e_value_type             = double;
    //!\brief The type of field::BITSCORE is fixed to double.
    using bitscore_type            = double;
    //!\brief The type of field::HEADER_PTR (default: alignment_file_header<typename traits_type::ref_ids>).
    using header_type              = alignment_file_header<typename traits_type::ref_ids>;

private:
    //!\brief The type of the aligned query sequence (second type of the pair of alignment_type).
    using alignment_query_type = std::conditional_t<
                                     selected_field_ids::contains(field::SEQ),
                                     gap_decorator<
                                         decltype(std::declval<sequence_type &>() | view::slice(0, 0))>,
                                     typename traits_type::template sequence_container<
                                         gapped<typename traits_type::sequence_alphabet>>>;

public:
    //!\brief The type of field::ALIGNMENT (default: std::pair<std::vector<gapped<dna5>>, std::vector<gapped<dna5>>>).
    using alignment_type = std::tuple<gap_decorator<ref_sequence_type>, alignment_query_type>;

    //!\brief The previously defined types aggregated in a seqan3::type_list.
    using field_types = type_list<sequence_type,
                                  id_type,
                                  offset_type,
                                  ref_sequence_type,
                                  ref_id_type,
                                  ref_offset_type,
                                  alignment_type,
                                  mapq_type,
                                  quality_type,
                                  flag_type,
                                  mate_type,
                                  sam_tag_dictionary,
                                  e_value_type,
                                  bitscore_type,
                                  header_type *>;

    /*!\brief The subset of seqan3::field tags that are valid for this file; order corresponds to the types in
     * \ref field_types.
     */
    using field_ids = fields<field::SEQ,
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
                             field::HEADER_PTR>;

    static_assert([] () constexpr
                  {
                      for (field f : selected_field_ids::as_array)
                          if (!field_ids::contains(f))
                              return false;
                      return true;
                  }(),
                  "You selected a field that is not valid for aligment files, please refer to the documentation "
                  "of alignment_file_input::field_ids for the accepted values.");

    //!\brief The type of the record, a specialisation of seqan3::record; acts as a tuple of the selected field types.
    using record_type = record<detail::select_types_with_ids_t<field_types, field_ids, selected_field_ids>,
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
    //!\brief The const_reference type is void because files are not const-iterable.
    using const_reference   = void;
    //!\brief An unsigned integer type, usually std::size_t.
    using size_type         = size_t;
    //!\brief A signed integer type, usually std::ptrdiff_t.
    using difference_type   = std::make_signed_t<size_t>;
    //!\brief The iterator type of this view (an input iterator).
    using iterator          = detail::in_file_iterator<alignment_file_input>;
    //!\brief The const iterator type is void because files are not const-iterable.
    using const_iterator    = void;
    //!\brief The type returned by end().
    using sentinel          = std::ranges::default_sentinel_t;
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    //!\brief Default constructor is explicitly deleted, you need to give a stream or file name.
    alignment_file_input() = delete;
    //!\brief Copy construction is explicitly deleted because you cannot have multiple access to the same file.
    alignment_file_input(alignment_file_input const &) = delete;
    //!\brief Copy assignment is explicitly deleted because you cannot have multiple access to the same file.
    alignment_file_input & operator=(alignment_file_input const &) = delete;
    //!\brief Move construction is defaulted.
    alignment_file_input(alignment_file_input &&) = default;
    //!\brief Move assignment is defaulted.
    alignment_file_input & operator=(alignment_file_input &&) = default;
    //!\brief Destructor is defaulted.
    ~alignment_file_input() = default;

    /*!\brief Construct from filename.
     * \param[in] filename    Path to the file you wish to open.
     * \param[in] fields_tag  A seqan3::fields tag. [optional]
     * \throws seqan3::file_open_error If the file could not be opened, e.g. non-existent, non-readable, unknown format.
     *
     * \details
     *
     * In addition to the file name, you may specify a custom seqan3::fields object
     * (e.g. `seqan3::fields<seqan3::field::SEQ>{}`) which may be easier than
     * defining all the template parameters.
     *
     * ### Decompression
     *
     * This constructor transparently applies a decompression stream on top of the file stream in case
     * the file is detected as being compressed.
     * See the section on \link io_compression compression and decompression \endlink for more information.
     */
    alignment_file_input(std::filesystem::path filename,
                         selected_field_ids const & SEQAN3_DOXYGEN_ONLY(fields_tag) = selected_field_ids{}) :
        primary_stream{new std::ifstream{filename, std::ios_base::in | std::ios::binary}, stream_deleter_default}
    {
        init(filename);
    }

    /*!\brief Construct from an existing stream and with specified format.
     * \tparam stream_t      The stream type; must model seqan3::IStream2.
     * \tparam file_format   The format of the file in the stream, must model seqan3::AlignmentFileInputFormat.
     * \param[in] stream     The stream to operate on; must be derived of std::basic_istream.
     * \param[in] format_tag The file format tag.
     * \param[in] fields_tag A seqan3::fields tag. [optional]
     *
     * \details
     *
     * In addition to the stream and the format, you may specify a custom seqan3::fields object
     * (e.g. `seqan3::fields<seqan3::field::SEQ>{}`) which may be easier than
     * defining all the template parameters.
     *
     * ### Decompression
     *
     * This constructor transparently applies a decompression stream on top of the stream in case
     * it is detected as being compressed.
     * See the section on \link io_compression compression and decompression \endlink for more information.
     */
    template <IStream2 stream_t, AlignmentFileInputFormat file_format>
    alignment_file_input(stream_t                 & stream,
                         file_format        const & SEQAN3_DOXYGEN_ONLY(format_tag),
                         selected_field_ids const & SEQAN3_DOXYGEN_ONLY(fields_tag) = selected_field_ids{}) :
        primary_stream{&stream, stream_deleter_noop}
    {
        init(file_format{});
    }

    //!\overload
    template <IStream2 stream_t, AlignmentFileInputFormat file_format>
    alignment_file_input(stream_t                && stream,
                         file_format        const & SEQAN3_DOXYGEN_ONLY(format_tag),
                         selected_field_ids const & SEQAN3_DOXYGEN_ONLY(fields_tag) = selected_field_ids{}) :
        primary_stream{new stream_t{std::move(stream)}, stream_deleter_default}
    {
        init(file_format{});
    }

    /*!\brief Construct from filename and given additional reference information.
     * \param[in] filename        Path to the file you wish to open.
     * \param[in] ref_ids         A range containing the reference ids that correspond to the SAM/BAM file.
     * \param[in] ref_sequences   A range containing the reference sequences that correspond to the SAM/BAM file.
     * \param[in] fields_tag      A seqan3::fields tag. [optional]
     * \throws seqan3::file_open_error If the file could not be opened, e.g. non-existent, non-readable, unknown format.
     *
     * \details
     *
     * The reference information given by the ids (names) and sequences will be used to construct a proper alignment
     * when reading in SAM or BAM files. If you are not interested in the full alignment, call the constructor without
     * the parameters.
     *
     * In addition to the file name and reference information, you may specify a custom seqan3::fields object
     * (e.g. `seqan3::fields<seqan3::field::SEQ>{}`) which may be easier than
     * defining all the template parameters.
     *
     * ### Decompression
     *
     * This constructor transparently applies a decompression stream on top of the file stream in case
     * the file is detected as being compressed.
     * See the section on \link io_compression compression and decompression \endlink for more information.
     */
    alignment_file_input(std::filesystem::path filename,
                         typename traits_type::ref_ids & ref_ids,
                         typename traits_type::ref_sequences & ref_sequences,
                         selected_field_ids const & SEQAN3_DOXYGEN_ONLY(fields_tag) = selected_field_ids{}) :
        primary_stream{new std::ifstream{filename, std::ios_base::in | std::ios::binary}, stream_deleter_default}
    {
        // initialize reference information
        set_references(ref_ids, ref_sequences);

        init(filename);
    }

    /*!\brief Construct from an existing stream and with specified format.
     * \tparam stream_t           The stream type; must model seqan3::IStream2.
     * \tparam    file_format     The format of the file in the stream; must model seqan3::AlignmentFileInputFormat.
     * \param[in] stream          The stream to operate on; must be derived of std::basic_istream.
     * \param[in] ref_ids         A range containing the reference ids that correspond to the SAM/BAM file.
     * \param[in] ref_sequences   A range containing the reference sequences that correspond to the SAM/BAM file.
     * \param[in] format_tag      The file format tag.
     * \param[in] fields_tag      A seqan3::fields tag. [optional]
     *
     * \details
     *
     * The reference information given by the ids (names) and sequences will be used to construct a proper alignment
     * when reading in SAM or BAM files. If you are not interested in the full alignment, you do not need to specify
     * those information.
     *
     * In addition to the stream, reference information and format, you may specify a custom seqan3::fields object
     * (e.g. `seqan3::fields<seqan3::field::SEQ>{}`) which may be easier than
     * defining all the template parameters.
     *
     * ### Decompression
     *
     * This constructor transparently applies a decompression stream on top of the stream in case
     * it is detected as being compressed.
     * See the section on \link io_compression compression and decompression \endlink for more information.
     */
    template <IStream2 stream_t, AlignmentFileInputFormat file_format>
    alignment_file_input(stream_t & stream,
                         typename traits_type::ref_ids & ref_ids,
                         typename traits_type::ref_sequences & ref_sequences,
                         file_format const & SEQAN3_DOXYGEN_ONLY(format_tag),
                         selected_field_ids const & SEQAN3_DOXYGEN_ONLY(fields_tag) = selected_field_ids{}) :
        primary_stream{&stream, stream_deleter_noop}
    {
        // initialize reference information
        set_references(ref_ids, ref_sequences);

        init(file_format{});
    }

    //!\overload
    template <IStream2 stream_t, AlignmentFileInputFormat file_format>
    alignment_file_input(stream_t && stream,
                         typename traits_type::ref_ids & ref_ids,
                         typename traits_type::ref_sequences & ref_sequences,
                         file_format const & SEQAN3_DOXYGEN_ONLY(format_tag),
                         selected_field_ids const & SEQAN3_DOXYGEN_ONLY(fields_tag) = selected_field_ids{}) :
        primary_stream{new stream_t{std::move(stream)}, stream_deleter_default}
    {
        // initialize reference information
        set_references(ref_ids, ref_sequences);

        init(file_format{});
    }

    //!\cond
    // explicitly delete rvalues for reference information
    alignment_file_input(std::filesystem::path,
                         typename traits_type::ref_ids &&,
                         typename traits_type::ref_sequences &&,
                         selected_field_ids const &) = delete;

    template <IStream2 stream_t, AlignmentFileInputFormat file_format>
    alignment_file_input(stream_t &&,
                         typename traits_type::ref_ids &&,
                         typename traits_type::ref_sequences &&,
                         file_format const &,
                         selected_field_ids const &) = delete;
    //!\endcond
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
     * \snippet test/snippet/io/alignment_file/alignment_file_input.cpp begin_and_front
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
     * and begin also always points to the current record on single pass input ranges:
     *
     * \snippet test/snippet/io/alignment_file/alignment_file_input.cpp begin_and_front
     *
     * In most situations using the iterator interface or a range-based for-loop are preferable to using front(),
     * because you can only move to the next record via the iterator.
     *
     * In any case, don't forget the reference! If you want to save the data from the record elsewhere, use move:
     *
     * \snippet test/snippet/io/alignment_file/alignment_file_input.cpp front
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

    //!\brief The options are public and its members can be set directly.
    alignment_file_input_options<typename traits_type::sequence_legal_alphabet> options;

    /*!\brief Access the file's header.
     *
     * \details
     *
     * You can access the header directly after the construction **with reference information** of the file object.
     *
     * ### Example
     *
     * \snippet test/snippet/io/alignment_file/alignment_file_input.cpp sam_file
     *
     * \snippet test/snippet/io/alignment_file/alignment_file_input.cpp get_header
     *
     * \sa seqan3::alignment_file_header
     */
    header_type & header()
    {
        return *header_ptr;
    }

protected:
    //!\privatesection

    //!/brief Initialisation based on a filename.
    void init(std::filesystem::path & filename)
    {
        // open stream
        if (!primary_stream->good())
            throw file_open_error{"Could not open file " + filename.string() + " for reading."};

        secondary_stream = detail::make_secondary_istream(*primary_stream, filename);
        detail::set_format(format, filename);

        // buffer first record
        read_next_record();
    }

    //!/brief Initialisation based on a format (construction via stream).
    template <typename format_type>
    void init(format_type const &)
    {
        static_assert(meta::in<valid_formats, format_type>::value,
                      "You selected a format that is not in the valid_formats of this file.");

        format = detail::alignment_file_input_format<format_type>{};
        secondary_stream = detail::make_secondary_istream(*primary_stream);

        // buffer first record
        read_next_record();
    }

    //!\brief The file header object.
    std::unique_ptr<header_type> header_ptr{new header_type{}};

    /*!\name Data buffers
     * \{
     */
    //!\brief Buffer for a single record.
    record_type record_buffer;
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

    //!\brief File is one position behind the last record.
    bool at_end{false};

    //!\brief Type of the format, an std::variant over the `valid_formats`.
    using format_type = typename detail::variant_from_tags<valid_formats, detail::alignment_file_input_format>::type;

    //!\brief The actual std::variant holding a pointer to the detected/selected format.
    format_type format;
    //!\}

    /*!\name Reference information
     * \{
     */
    //!\brief A pointer to the reference sequence information if given on construction.
    typename traits_type::ref_sequences const * reference_sequences_ptr{nullptr};

    /*!\brief Updates the reference information members and the header.
     *
     * \details
     *
     * The SAM format only provides semi-alignments because the reference sequence
     * is not stored explicitly. In order to be able to read in full alignments,
     * additional reference information can be given to the alignment file on construction.
     * Note that the reference ids (names) must correspond to the exact spelling
     * in the SAM/BAM file otherwise an exception will be thrown when reading.
     */
    template <std::ranges::ForwardRange ref_sequences_t>
    void set_references(typename traits_type::ref_ids & ref_ids, ref_sequences_t && ref_sequences)
    {
        assert(std::ranges::distance(ref_ids) == std::ranges::distance(ref_sequences));

        header_ptr = std::unique_ptr<header_type>{std::make_unique<header_type>(ref_ids)};
        reference_sequences_ptr = &ref_sequences;

        // initialise reference map and ref_dict if ref_ids are non-empty
        for (int32_t idx = 0; idx < std::ranges::distance(ref_ids); ++idx)
        {
            header_ptr->ref_id_info.emplace_back(std::ranges::distance(ref_sequences[idx]), "");

            if constexpr (std::ranges::ContiguousRange<reference_t<typename traits_type::ref_ids>> &&
                          std::ranges::SizedRange<reference_t<typename traits_type::ref_ids>> &&
                          ForwardingRange<reference_t<typename traits_type::ref_ids>>)
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
    //!\}

    //!\brief Tell the format to move to the next record and update the buffer.
    void read_next_record()
    {
        // clear the record
        record_buffer.clear();
        detail::get_or_ignore<field::HEADER_PTR>(record_buffer) = header_ptr.get();

        // at end if we could not read further
        if (std::istreambuf_iterator<stream_char_type>{*secondary_stream} ==
            std::istreambuf_iterator<stream_char_type>{})
        {
            at_end = true;
            return;
        }

        auto call_read_func = [this] (auto & ref_seq_info)
        {
            std::visit([&] (auto & f)
            {
                f.read(*secondary_stream,
                       options,
                       ref_seq_info,
                       *header_ptr,
                       detail::get_or_ignore<field::SEQ>(record_buffer),
                       detail::get_or_ignore<field::QUAL>(record_buffer),
                       detail::get_or_ignore<field::ID>(record_buffer),
                       detail::get_or_ignore<field::OFFSET>(record_buffer),
                       detail::get_or_ignore<field::REF_SEQ>(record_buffer),
                       detail::get_or_ignore<field::REF_ID>(record_buffer),
                       detail::get_or_ignore<field::REF_OFFSET>(record_buffer),
                       detail::get_or_ignore<field::ALIGNMENT>(record_buffer),
                       detail::get_or_ignore<field::FLAG>(record_buffer),
                       detail::get_or_ignore<field::MAPQ>(record_buffer),
                       detail::get_or_ignore<field::MATE>(record_buffer),
                       detail::get_or_ignore<field::TAGS>(record_buffer),
                       detail::get_or_ignore<field::EVALUE>(record_buffer),
                       detail::get_or_ignore<field::BIT_SCORE>(record_buffer));

            }, format);
        };

        assert(!format.valueless_by_exception());

        if constexpr (!std::Same<typename traits_type::ref_sequences, ref_info_not_given>)
            call_read_func(*reference_sequences_ptr);
        else
            call_read_func(std::ignore);
    }

    //!\brief Befriend iterator so it can access the buffers.
    friend iterator;
};

/*!\name Type deduction guides
 * \relates seqan3::alignment_file_input
 * \{
 */
//!\brief Deduce selected fields, file_format and stream char type, default the rest.
template <IStream2                 stream_type,
          AlignmentFileInputFormat file_format,
          detail::Fields           selected_field_ids>
alignment_file_input(stream_type && stream,
                     file_format const &,
                     selected_field_ids const &)
    -> alignment_file_input<typename alignment_file_input<>::traits_type,       // actually use the default
                            selected_field_ids,
                            type_list<file_format>,
                            typename std::remove_reference_t<stream_type>::char_type>;

//!\brief Deduce selected fields, file_format and stream char type, default the rest.
template <IStream2                 stream_type,
          AlignmentFileInputFormat file_format,
          detail::Fields           selected_field_ids>
alignment_file_input(stream_type & stream,
                     file_format const &,
                     selected_field_ids const &)
    -> alignment_file_input<typename alignment_file_input<>::traits_type,       // actually use the default
                            selected_field_ids,
                            type_list<file_format>,
                            typename std::remove_reference_t<stream_type>::char_type>;

//!\brief Deduce file_format and stream char type, default the rest.
template <IStream2                 stream_type,
          AlignmentFileInputFormat file_format>
alignment_file_input(stream_type && stream,
                     file_format const &)
    -> alignment_file_input<typename alignment_file_input<>::traits_type,        // actually use the default
                            typename alignment_file_input<>::selected_field_ids, // actually use the default
                            type_list<file_format>,
                            typename std::remove_reference_t<stream_type>::char_type>;

//!\brief Deduce file_format and stream char type, default the rest.
template <IStream2                 stream_type,
          AlignmentFileInputFormat file_format>
alignment_file_input(stream_type & stream,
                     file_format const &)
    -> alignment_file_input<typename alignment_file_input<>::traits_type,        // actually use the default
                            typename alignment_file_input<>::selected_field_ids, // actually use the default
                            type_list<file_format>,
                            typename std::remove_reference_t<stream_type>::char_type>;

//!\brief Deduce selected fields, ref_sequences_t and ref_ids_t, default the rest.
template <std::ranges::ForwardRange           ref_ids_t,
          std::ranges::ForwardRange           ref_sequences_t,
          detail::Fields                      selected_field_ids>
alignment_file_input(std::filesystem::path path,
                     ref_ids_t &,
                     ref_sequences_t &,
                     selected_field_ids const &)
    -> alignment_file_input<alignment_file_input_default_traits<std::remove_reference_t<ref_sequences_t>,
                                                                std::remove_reference_t<ref_ids_t>>,
                            selected_field_ids,
                            typename alignment_file_input<>::valid_formats,      // actually use the default
                            typename alignment_file_input<>::stream_char_type>;  // actually use the default

//!\brief Deduce ref_sequences_t and ref_ids_t, default the rest.
template <std::ranges::ForwardRange ref_ids_t,
          std::ranges::ForwardRange ref_sequences_t>
alignment_file_input(std::filesystem::path path,
                     ref_ids_t &,
                     ref_sequences_t &)
    -> alignment_file_input<alignment_file_input_default_traits<std::remove_reference_t<ref_sequences_t>,
                                                                std::remove_reference_t<ref_ids_t>>,
                            typename alignment_file_input<>::selected_field_ids, // actually use the default
                            typename alignment_file_input<>::valid_formats,      // actually use the default
                            typename alignment_file_input<>::stream_char_type>;  // actually use the default

//!\brief Deduce selected fields, ref_sequences_t and ref_ids_t, file format and stream char type.
template <IStream2                  stream_type,
          std::ranges::ForwardRange ref_ids_t,
          std::ranges::ForwardRange ref_sequences_t,
          AlignmentFileInputFormat  file_format,
          detail::Fields            selected_field_ids>
alignment_file_input(stream_type && stream,
                     ref_ids_t &,
                     ref_sequences_t &,
                     file_format const &,
                     selected_field_ids const &)
    -> alignment_file_input<alignment_file_input_default_traits<std::remove_reference_t<ref_sequences_t>,
                                                                std::remove_reference_t<ref_ids_t>>,
                            selected_field_ids,
                            type_list<file_format>,
                            typename std::remove_reference_t<stream_type>::char_type>;

//!\brief Deduce selected fields, ref_sequences_t and ref_ids_t, file format and stream char type.
template <IStream2                  stream_type,
          std::ranges::ForwardRange ref_ids_t,
          std::ranges::ForwardRange ref_sequences_t,
          AlignmentFileInputFormat  file_format,
          detail::Fields            selected_field_ids>
alignment_file_input(stream_type & stream,
                     ref_ids_t &,
                     ref_sequences_t &,
                     file_format const &,
                     selected_field_ids const &)
    -> alignment_file_input<alignment_file_input_default_traits<std::remove_reference_t<ref_sequences_t>,
                                                                std::remove_reference_t<ref_ids_t>>,
                            selected_field_ids,
                            type_list<file_format>,
                            typename std::remove_reference_t<stream_type>::char_type>;

//!\brief Deduce ref_sequences_t and ref_ids_t, file format and stream char type.
template <IStream2                  stream_type,
          std::ranges::ForwardRange ref_ids_t,
          std::ranges::ForwardRange ref_sequences_t,
          AlignmentFileInputFormat  file_format>
alignment_file_input(stream_type && stream,
                     ref_ids_t &,
                     ref_sequences_t &,
                     file_format const &)
    -> alignment_file_input<alignment_file_input_default_traits<std::remove_reference_t<ref_sequences_t>,
                                                                std::remove_reference_t<ref_ids_t>>,
                            typename alignment_file_input<>::selected_field_ids, // actually use the default
                            type_list<file_format>,
                            typename std::remove_reference_t<stream_type>::char_type>;

//!\brief Deduce selected fields, ref_sequences_t and ref_ids_t, file format and stream char type.
template <IStream2                  stream_type,
          std::ranges::ForwardRange ref_ids_t,
          std::ranges::ForwardRange ref_sequences_t,
          AlignmentFileInputFormat  file_format>
alignment_file_input(stream_type & stream,
                     ref_ids_t &,
                     ref_sequences_t &,
                     file_format const &)
    -> alignment_file_input<alignment_file_input_default_traits<std::remove_reference_t<ref_sequences_t>,
                                                                std::remove_reference_t<ref_ids_t>>,
                            typename alignment_file_input<>::selected_field_ids, // actually use the default
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
 * \ingroup alignment_file
 * \see std::tuple_size_v
 */
template <seqan3::AlignmentFileInputTraits                    traits_type,
          seqan3::detail::Fields                              selected_field_ids,
          seqan3::detail::TypeListOfAlignmentFileInputFormats valid_formats,
          std::Integral                                       stream_char_t>
struct tuple_size<seqan3::alignment_file_input<traits_type, selected_field_ids, valid_formats, stream_char_t>>
{
    //!\brief The value equals the number of selected fields in the file.
    static constexpr size_t value = selected_field_ids::as_array.size();
};

/*!\brief Obtains the type of the specified element.
 * \implements seqan3::TransformationTrait
 * \ingroup alignment_file
 * \see [std::tuple_element](https://en.cppreference.com/w/cpp/utility/tuple/tuple_element)
 */
template <size_t                                              elem_no,
          seqan3::AlignmentFileInputTraits                    traits_type,
          seqan3::detail::Fields                              selected_field_ids,
          seqan3::detail::TypeListOfAlignmentFileInputFormats valid_formats,
          std::Integral                                       stream_char_t>
struct tuple_element<elem_no, seqan3::alignment_file_input<traits_type, selected_field_ids, valid_formats, stream_char_t>>
    : tuple_element<elem_no, typename seqan3::alignment_file_input<traits_type,
                                                               selected_field_ids,
                                                               valid_formats,
                                                               stream_char_t>::file_as_tuple_type>
{};

} // namespace std

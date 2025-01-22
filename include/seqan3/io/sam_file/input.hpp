// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::sam_file_input and corresponding traits classes.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <cassert>
#include <concepts>
#include <filesystem>
#include <fstream>
#include <ranges>
#include <string>
#include <variant>
#include <vector>

#include <seqan3/alphabet/adaptation/char.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/cigar/cigar.hpp>
#include <seqan3/alphabet/nucleotide/dna15.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/alphabet/quality/qualified.hpp>
#include <seqan3/io/detail/in_file_iterator.hpp>
#include <seqan3/io/detail/misc_input.hpp>
#include <seqan3/io/detail/record.hpp>
#include <seqan3/io/exception.hpp>
#include <seqan3/io/sam_file/format_bam.hpp>
#include <seqan3/io/sam_file/format_sam.hpp>
#include <seqan3/io/sam_file/input_format_concept.hpp>
#include <seqan3/io/sam_file/record.hpp>
#include <seqan3/io/sam_file/sam_flag.hpp>
#include <seqan3/io/stream/concept.hpp>
#include <seqan3/utility/tuple/concept.hpp>
#include <seqan3/utility/type_list/traits.hpp>
#include <seqan3/utility/type_traits/detail/transformation_trait_or.hpp>
#include <seqan3/utility/views/repeat_n.hpp>
#include <seqan3/utility/views/slice.hpp>

namespace seqan3
{

// ---------------------------------------------------------------------------------------------------------------------
// sam_file_input_traits
// ---------------------------------------------------------------------------------------------------------------------

/*!\interface seqan3::sam_file_input_traits <>
 * \brief The requirements a traits_type for seqan3::sam_file_input must meet.
 * \ingroup io_sam_file
 */
/*!\name Requirements for seqan3::sam_file_input_traits
 * \brief You can expect these **member types** of all types that model seqan3::sam_file_input_traits.
 * \memberof seqan3::sam_file_input_traits
 * \{
 */
/*!\typedef using sequence_alphabet
 * \brief Alphabet of the characters for the seqan3::field::seq; must model seqan3::alphabet.
 */
/*!\typedef using sequence_legal_alphabet
 * \brief Intermediate alphabet for seqan3::field::seq; must model seqan3::alphabet and be convertible to
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
 * must model seqan3::sequence_container.
 */
/*!\typedef using id_container
 * \brief Type template of the seqan3::field::id, a container template over `char`;
 * must model seqan3::sequence_container.
 */
/*!\typedef using quality_alphabet
 * \brief Alphabet of the characters for the seqan3::field::qual; must model seqan3::writable_quality_alphabet.
 */
/*!\typedef using quality_container
 * \brief Type template of the seqan3::field::qual, a container template over `quality_alphabet`;
 * must model seqan3::sequence_container.
 */
/*!\typedef using ref_sequences
 * \brief The type of range over reference sequences; must model std::ranges::forward_range,
 *        the value_type must also model std::ranges::forward_range, and the value type of the value type
 *        must model seqan3::alphabet (e.g. std::vector<std::vector<dna4>>).
 *
 * \attention This type is the first template parameter of seqan3::sam_file_input_default_traits and should not be
 *            manually configured in order to allow for automatic type deduction from reference information input on
 *            construction.
 */
/*!\typedef using ref_ids
 * \brief The type of range over reference sequences; must model std::ranges::forward_range,
 *        the value_type must also model std::ranges::forward_range, and the value type of the value type
 *        must model seqan3::alphabet (e.g. std::vector<string>).
 *
 * \attention This type is the second template parameter of seqan3::sam_file_input_default_traits and should not be
 *            manually configured in order to allow for automatic type deduction from reference information input on
 *            construction.
 */
//!\}
//!\cond
template <typename t>
concept sam_file_input_traits = requires (t v) {
    // field::seq
    requires writable_alphabet<typename t::sequence_alphabet>;
    requires writable_alphabet<typename t::sequence_legal_alphabet>;
    requires explicitly_convertible_to<typename t::sequence_legal_alphabet, typename t::sequence_alphabet>;
    requires sequence_container<typename t::template sequence_container<typename t::sequence_alphabet>>;

    // field::id
    requires sequence_container<typename t::template id_container<char>>;

    // field::qual
    requires writable_quality_alphabet<typename t::quality_alphabet>;
    requires sequence_container<typename t::template quality_container<typename t::quality_alphabet>>;

    // field::ref_seq
    // either ref_info_not_given or a range over ranges over alphabet (e.g. std::vector<dna4_vector>)
    requires std::same_as<typename t::ref_sequences, ref_info_not_given> || requires () {
        requires alphabet<std::ranges::range_reference_t<std::ranges::range_reference_t<typename t::ref_sequences>>>;
    };

    // field::ref_id
    requires alphabet<std::ranges::range_reference_t<std::ranges::range_reference_t<typename t::ref_ids>>>
                 && (!std::same_as<typename t::ref_sequences, ref_info_not_given>
                     || writable_alphabet<
                         std::ranges::range_reference_t<std::ranges::range_reference_t<typename t::ref_ids>>>);
    requires std::ranges::forward_range<std::ranges::range_reference_t<typename t::ref_ids>>;
    requires std::ranges::forward_range<typename t::ref_ids>;

    // field::ref_offset is fixed to std::optional<int32_t>
    // field::flag is fixed to seqan3::sam_flag
    // field::mapq is fixed to uint8_t
    // field::evalue is fixed to double
    // field::bitscore is fixed to double
    // field::mate is fixed to std::tuple<ref_id_container<ref_id_alphabet>, ref_offset_type, int32_t>
};
//!\endcond

// ---------------------------------------------------------------------------------------------------------------------
// sam_file_input_default_traits
// ---------------------------------------------------------------------------------------------------------------------

/*!\brief The default traits for seqan3::sam_file_input
 * \implements sam_file_input_traits
 * \ingroup io_sam_file
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
 * \include test/snippet/io/sam_file/sam_file_input_my_traits.cpp
 */
template <typename ref_sequences_t = ref_info_not_given, typename ref_ids_t = std::deque<std::string>>
struct sam_file_input_default_traits
{
    /*!\name Member types
     * \brief Definitions to model seqan3::sam_file_input_traits.
     * \{
     */

    //!\brief The sequence alphabet is seqan3::dna5.
    using sequence_alphabet = dna5;

    //!\brief The legal sequence alphabet for parsing is seqan3::dna15.
    using sequence_legal_alphabet = dna15;

    //!\brief The container for a sequence is std::vector.
    template <typename _sequence_alphabet>
    using sequence_container = std::vector<_sequence_alphabet>;

    //!\brief The string type for an identifier is std::basic_string.
    template <typename _id_alphabet>
    using id_container = std::basic_string<_id_alphabet>;

    //!\brief The alphabet for a quality annotation is seqan3::phred42.
    using quality_alphabet = phred42;

    //!\brief The string type for a quality annotation is std::vector.
    template <typename _quality_alphabet>
    using quality_container = std::vector<_quality_alphabet>;

    //!\brief The type of the reference sequences is deduced on construction.
    using ref_sequences = ref_sequences_t;

    //!\brief The type of the reference identifiers is deduced on construction.
    using ref_ids = ref_ids_t;
    //!\}
};

// ---------------------------------------------------------------------------------------------------------------------
// sam_file_input
// ---------------------------------------------------------------------------------------------------------------------

/*!\brief A class for reading SAM files, both SAM and its binary representation BAM are supported.
 * \ingroup io_sam_file
 * \tparam traits_type          An auxiliary type that defines certain member types and constants, must model
 *                              seqan3::sam_file_input_traits.
 * \tparam selected_field_ids   A seqan3::fields type with the list and order of desired record entries; all fields
 *                              must be in seqan3::sam_file_input::field_ids.
 * \tparam valid_formats        A seqan3::type_list of the selectable formats (each must meet
 *                              seqan3::sam_file_input_format).
 *
 * \details
 *
 * \include{doc} doc/fragments/io_sam_file_input.md
 *
 * \remark For a complete overview, take a look at \ref io_sam_file
 */
template <sam_file_input_traits traits_type_ = sam_file_input_default_traits<>,
          detail::fields_specialisation selected_field_ids_ = fields<field::seq,
                                                                     field::id,
                                                                     field::ref_id,
                                                                     field::ref_offset,
                                                                     field::cigar,
                                                                     field::mapq,
                                                                     field::qual,
                                                                     field::flag,
                                                                     field::mate,
                                                                     field::tags,
                                                                     field::header_ptr>,
          detail::type_list_of_sam_file_input_formats valid_formats_ = type_list<format_sam, format_bam>>
class sam_file_input
{
public:
    /*!\name Template arguments
     * \brief Exposed as member types for public access.
     * \{
     */
    //!\brief A traits type that defines aliases and template for storage of the fields.
    using traits_type = traits_type_;
    //!\brief A seqan3::fields list with the fields selected for the record.
    using selected_field_ids = selected_field_ids_;
    //!\brief A seqan3::type_list with the possible formats.
    using valid_formats = valid_formats_;
    //!\brief Character type of the stream(s).
    using stream_char_type = char;
    //!\}

private:
    //!\brief The dummy ref sequence type if no reference information were given.
    using dummy_ref_type = decltype(views::repeat_n(typename traits_type::sequence_alphabet{}, size_t{})
                                    | std::views::transform(detail::access_restrictor_fn{}));

    //!\brief The unsliced ref sequence type if reference information were given.
    using ref_sequence_unsliced_type = detail::lazy_conditional_t<
        std::ranges::range<typename traits_type::ref_sequences const>,
        detail::lazy<std::ranges::range_reference_t, typename traits_type::ref_sequences const>,
        dummy_ref_type>;

    //!\brief The ref sequence type if reference information were given.
    using ref_sequence_sliced_type = decltype(std::declval<ref_sequence_unsliced_type>() | views::slice(0, 0));

public:
    /*!\name Field types and record type
     * \brief These types are relevant for record/row-based reading; they may be manipulated via the \ref traits_type
     * to achieve different storage behaviour.
     * \{
     */
    //!\brief The type of field::seq (default std::vector<seqan3::dna5>).
    using sequence_type = typename traits_type::template sequence_container<typename traits_type::sequence_alphabet>;
    //!\brief The type of field::id (default std::string by default).
    using id_type = typename traits_type::template id_container<char>;
    /*!\brief The type of field::ref_seq (default depends on construction).
     *
     * If no reference information are given on construction, this type deduces to a sized view that throws on
     * access (since there is nothing to access anyway). If the reference information are given, the type is deduced
     * to a view over the given input reference sequence type such that no sequence information is copied.
     */
    using ref_sequence_type = std::conditional_t<std::same_as<typename traits_type::ref_sequences, ref_info_not_given>,
                                                 dummy_ref_type,
                                                 ref_sequence_sliced_type>;
    /*!\brief The type of field::ref_id is fixed to std::optional<int32_t>.
     *
     * To be consistent with the BAM format, the field::ref_id will hold the index to the actual reference
     * information stored in the header. If a read is unmapped, the optional will remain valueless.
     *
     * \attention SeqaAn3 transforms the 1-based SAM format position into a 0-based position.
     */
    using ref_id_type = std::optional<int32_t>;
    /*!\brief The type of field::ref_offset is fixed to a std::optional<int32_t>.
     *
     * The SAM format is 1-based and a 0 in the ref_offset field indicated an unmapped read. Since we convert 1-based
     * positions to 0-based positions when reading the SAM format, we model the ref_offset_type as a std::optional.
     * If the input value is 0, the std::optional will remain valueless.
     */
    using ref_offset_type = std::optional<int32_t>;
    //!\brief The type of field::mapq is fixed to uint8_t.
    using mapq_type = uint8_t;
    //!\brief The type of field::qual (default std::vector<seqan3::phred42>).
    using quality_type = typename traits_type::template quality_container<typename traits_type::quality_alphabet>;
    //!\brief The type of field::flag is fixed to seqan3::sam_flag.
    using flag_type = sam_flag;
    //!\brief The type of field::cigar is fixed to std::vector<cigar>.
    using cigar_type = std::vector<cigar>;
    //!\brief The type of field::mate is fixed to std::tuple<ref_id_type, ref_offset_type, int32_t>).
    using mate_type = std::tuple<ref_id_type, ref_offset_type, int32_t>;
    //!\brief The type of field::header_ptr (default: sam_file_header<typename traits_type::ref_ids>).
    using header_type = sam_file_header<typename traits_type::ref_ids>;

    //!\brief The previously defined types aggregated in a seqan3::type_list.
    using field_types = type_list<sequence_type,
                                  id_type,
                                  ref_id_type,
                                  ref_offset_type,
                                  std::vector<cigar>,
                                  mapq_type,
                                  quality_type,
                                  flag_type,
                                  mate_type,
                                  sam_tag_dictionary,
                                  header_type *>;

    /*!\brief The subset of seqan3::field tags valid for this file; order corresponds to the types in \ref field_types.
     *
     * The SAM file abstraction supports reading 10 different fields:
     *
     *   1. seqan3::field::seq
     *   2. seqan3::field::id
     *   3. seqan3::field::ref_id
     *   4. seqan3::field::ref_offset
     *   5. seqan3::field::cigar
     *   6. seqan3::field::mapq
     *   7. seqan3::field::qual
     *   8. seqan3::field::flag
     *   9. seqan3::field::mate
     *   10. seqan3::field::tags
     *
     * There exists one more field for SAM files, the seqan3::field::header_ptr, but this field is mostly used
     * internally. Please see the seqan3::sam_file_output::header member function for details on how to access
     * the seqan3::sam_file_header of the file.
     */
    using field_ids = fields<field::seq,
                             field::id,
                             field::ref_id,
                             field::ref_offset,
                             field::cigar,
                             field::mapq,
                             field::qual,
                             field::flag,
                             field::mate,
                             field::tags,
                             field::header_ptr>;

    static_assert(!selected_field_ids::contains(field::alignment),
                  "The seqan3::field::alignment was removed from the allowed fields for seqan3::sam_file_input. "
                  "Only seqan3::field::cigar is supported. Please see seqan3::alignment_from_cigar on how to get an "
                  "alignment from the cigar information.");

    static_assert(!selected_field_ids::contains(field::offset),
                  "The field::offset is deprecated. Please access field::cigar and retrieve the soft clipping (S) "
                  "value at the front of the CIGAR string (offset = 0 if there is no soft clipping at the front).");

    static_assert(
        []() constexpr
        {
            for (field f : selected_field_ids::as_array)
                if (!field_ids::contains(f))
                    return false;
            return true;
        }(),
        "You selected a field that is not valid for SAM files, please refer to the documentation "
        "of sam_file_input::field_ids for the accepted values.");

    //!\brief The type of the record, a specialisation of seqan3::record; acts as a tuple of the selected field types.
    using record_type =
        sam_record<detail::select_types_with_ids_t<field_types, field_ids, selected_field_ids>, selected_field_ids>;
    //!\}

    /*!\name Range associated types
     * \brief The types necessary to facilitate the behaviour of an input range (used in record-wise reading).
     * \{
     */
    //!\brief The value_type is the \ref record_type.
    using value_type = record_type;
    //!\brief The reference type.
    using reference = record_type &;
    //!\brief The const_reference type is void because files are not const-iterable.
    using const_reference = void;
    //!\brief An unsigned integer type, usually std::size_t.
    using size_type = size_t;
    //!\brief A signed integer type, usually std::ptrdiff_t.
    using difference_type = std::make_signed_t<size_t>;
    //!\brief The iterator type of this view (an input iterator).
    using iterator = detail::in_file_iterator<sam_file_input>;
    //!\brief The const iterator type is void because files are not const-iterable.
    using const_iterator = void;
    //!\brief The type returned by end().
    using sentinel = std::default_sentinel_t;
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    //!\brief Default constructor is explicitly deleted, you need to give a stream or file name.
    sam_file_input() = delete;
    //!\brief Copy construction is explicitly deleted because you cannot have multiple access to the same file.
    sam_file_input(sam_file_input const &) = delete;
    //!\brief Copy assignment is explicitly deleted because you cannot have multiple access to the same file.
    sam_file_input & operator=(sam_file_input const &) = delete;
    //!\brief Move construction is defaulted.
    sam_file_input(sam_file_input &&) = default;
    //!\brief Move assignment is defaulted.
    sam_file_input & operator=(sam_file_input &&) = default;
    //!\brief Destructor is defaulted.
    ~sam_file_input() = default;

    /*!\brief Construct from filename.
     * \param[in] filename    Path to the file you wish to open.
     * \param[in] fields_tag  A seqan3::fields tag. [optional]
     * \throws seqan3::file_open_error If the file could not be opened, e.g. non-existent, non-readable, unknown format.
     *
     * \details
     *
     * In addition to the file name, you may specify a custom seqan3::fields object
     * (e.g. `seqan3::fields<seqan3::field::seq>{}`) which may be easier than
     * defining all the template parameters.
     *
     * ### Decompression
     *
     * This constructor transparently applies a decompression stream on top of the file stream in case
     * the file is detected as being compressed.
     * See the section on \link io_compression compression and decompression \endlink for more information.
     */
    sam_file_input(std::filesystem::path filename,
                   selected_field_ids const & SEQAN3_DOXYGEN_ONLY(fields_tag) = selected_field_ids{}) :
        primary_stream{new std::ifstream{}, stream_deleter_default}
    {
        init_by_filename(std::move(filename));
    }

    /*!\brief Construct from an existing stream and with specified format.
     * \tparam stream_t      The stream type; must model seqan3::input_stream.
     * \tparam file_format   The format of the file in the stream, must model seqan3::sam_file_input_format.
     * \param[in] stream     The stream to operate on; must be derived of std::basic_istream.
     * \param[in] format_tag The file format tag.
     * \param[in] fields_tag A seqan3::fields tag. [optional]
     *
     * \details
     *
     * In addition to the stream and the format, you may specify a custom seqan3::fields object
     * (e.g. `seqan3::fields<seqan3::field::seq>{}`) which may be easier than
     * defining all the template parameters.
     *
     * ### Decompression
     *
     * This constructor transparently applies a decompression stream on top of the stream in case
     * it is detected as being compressed.
     * See the section on \link io_compression compression and decompression \endlink for more information.
     */
    template <input_stream stream_t, sam_file_input_format file_format>
        requires std::same_as<typename std::remove_reference_t<stream_t>::char_type, stream_char_type>
    sam_file_input(stream_t & stream,
                   file_format const & SEQAN3_DOXYGEN_ONLY(format_tag),
                   selected_field_ids const & SEQAN3_DOXYGEN_ONLY(fields_tag) = selected_field_ids{}) :
        primary_stream{&stream, stream_deleter_noop}
    {
        init_by_format<file_format>();
    }

    //!\overload
    template <input_stream stream_t, sam_file_input_format file_format>
        requires std::same_as<typename std::remove_reference_t<stream_t>::char_type, stream_char_type>
    sam_file_input(stream_t && stream,
                   file_format const & SEQAN3_DOXYGEN_ONLY(format_tag),
                   selected_field_ids const & SEQAN3_DOXYGEN_ONLY(fields_tag) = selected_field_ids{}) :
        primary_stream{new stream_t{std::move(stream)}, stream_deleter_default}
    {
        init_by_format<file_format>();
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
     * ### Reference information
     *
     * The reference information given by the IDs (names) and sequences will be used to keep the record entry
     * `seqan3::sam_file_input::record_type::reference_id()` consistent with the order imposed by `ref_ids`.
     * This way, you can use the value of `seqan3::sam_file_input::record_type::reference_id()` to access the lists
     * `ref_ids` and `ref_sequences` to retrieve the correct information for the current record.
     *
     * ### Selecting custom fields
     *
     * In addition to the file name and reference information, you may specify a custom seqan3::fields object
     * (e.g. `seqan3::fields<seqan3::field::seq>{}`) which may be easier than
     * defining all the template parameters.
     *
     * ### Decompression
     *
     * This constructor transparently applies a decompression stream on top of the file stream in case
     * the file is detected as being compressed.
     * See the section on \link io_compression compression and decompression \endlink for more information.
     */
    sam_file_input(std::filesystem::path filename,
                   typename traits_type::ref_ids & ref_ids,
                   typename traits_type::ref_sequences & ref_sequences,
                   selected_field_ids const & SEQAN3_DOXYGEN_ONLY(fields_tag) = selected_field_ids{}) :
        primary_stream{new std::ifstream{}, stream_deleter_default}
    {
        // initialize reference information
        set_references(ref_ids, ref_sequences);

        init_by_filename(std::move(filename));
    }

    /*!\brief Construct from an existing stream and with specified format.
     * \tparam stream_t           The stream type; must model seqan3::input_stream.
     * \tparam    file_format     The format of the file in the stream; must model seqan3::sam_file_input_format.
     * \param[in] stream          The stream to operate on; must be derived of std::basic_istream.
     * \param[in] ref_ids         A range containing the reference ids that correspond to the SAM/BAM file.
     * \param[in] ref_sequences   A range containing the reference sequences that correspond to the SAM/BAM file.
     * \param[in] format_tag      The file format tag.
     * \param[in] fields_tag      A seqan3::fields tag. [optional]
     *
     * \details
     *
     * ### Reference information
     *
     * The reference information given by the IDs (names) and sequences will be used to keep the record entry
     * `seqan3::sam_file_input::record_type::reference_id()` consistent with the order imposed by `ref_ids`.
     * This way, you can use the value of `seqan3::sam_file_input::record_type::reference_id()` to access the lists
     * `ref_ids` and `ref_sequences` to retrieve the correct information for the current record.
     *
     * ### Selecting custom fields
     *
     * In addition to the stream, reference information and format, you may specify a custom seqan3::fields object
     * (e.g. `seqan3::fields<seqan3::field::seq>{}`) which may be easier than
     * defining all the template parameters.
     *
     * ### Decompression
     *
     * This constructor transparently applies a decompression stream on top of the stream in case
     * it is detected as being compressed.
     * See the section on \link io_compression compression and decompression \endlink for more information.
     */
    template <input_stream stream_t, sam_file_input_format file_format>
    sam_file_input(stream_t & stream,
                   typename traits_type::ref_ids & ref_ids,
                   typename traits_type::ref_sequences & ref_sequences,
                   file_format const & SEQAN3_DOXYGEN_ONLY(format_tag),
                   selected_field_ids const & SEQAN3_DOXYGEN_ONLY(fields_tag) = selected_field_ids{}) :
        primary_stream{&stream, stream_deleter_noop}
    {
        // initialize reference information
        set_references(ref_ids, ref_sequences);

        init_by_format<file_format>();
    }

    //!\overload
    template <input_stream stream_t, sam_file_input_format file_format>
    sam_file_input(stream_t && stream,
                   typename traits_type::ref_ids & ref_ids,
                   typename traits_type::ref_sequences & ref_sequences,
                   file_format const & SEQAN3_DOXYGEN_ONLY(format_tag),
                   selected_field_ids const & SEQAN3_DOXYGEN_ONLY(fields_tag) = selected_field_ids{}) :
        primary_stream{new stream_t{std::move(stream)}, stream_deleter_default}
    {
        // initialize reference information
        set_references(ref_ids, ref_sequences);

        init_by_format<file_format>();
    }

    //!\cond
    // explicitly delete rvalues for reference information
    sam_file_input(std::filesystem::path,
                   typename traits_type::ref_ids &&,
                   typename traits_type::ref_sequences &&,
                   selected_field_ids const &) = delete;

    template <input_stream stream_t, sam_file_input_format file_format>
    sam_file_input(stream_t &&,
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
     * \throws seqan3::format_error
     *
     * Equals end() if the file is at end.
     *
     * \include test/snippet/io/sam_file/sam_file_input_begin_and_front.cpp
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
     * and begin also always points to the current record on single pass input ranges:
     *
     * \include test/snippet/io/sam_file/sam_file_input_begin_and_front.cpp
     *
     * In most situations using the iterator interface or a range-based for-loop are preferable to using front(),
     * because you can only move to the next record via the iterator.
     *
     * In any case, don't forget the reference! If you want to save the data from the record elsewhere, use move:
     *
     * \include test/snippet/io/sam_file/sam_file_input_front.cpp
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

    //!\brief The options are public and its members can be set directly.
    sam_file_input_options<typename traits_type::sequence_legal_alphabet> options;

    /*!\brief Access the file's header.
     *
     * \details
     *
     * You can access the header directly after the construction **with reference information** of the file object.
     *
     * ### Example
     *
     * \include test/snippet/io/sam_file/sam_file_input_get_header.cpp
     *
     * \sa seqan3::sam_file_header
     */
    header_type & header()
    {
        // make sure header is read
        if (!first_record_was_read)
        {
            read_next_record();
            first_record_was_read = true;
        }

        return *header_ptr;
    }

protected:
    //!\privatesection

    //!/brief Initialisation based on a filename.
    void init_by_filename(std::filesystem::path filename)
    {
        primary_stream->rdbuf()->pubsetbuf(stream_buffer.data(), stream_buffer.size());
        static_cast<std::basic_ifstream<char> *>(primary_stream.get())
            ->open(filename, std::ios_base::in | std::ios::binary);
        // open stream
        if (!primary_stream->good())
            throw file_open_error{"Could not open file " + filename.string() + " for reading."};

        secondary_stream = detail::make_secondary_istream(*primary_stream, filename);
        detail::set_format(format, filename);
    }

    //!/brief Initialisation based on a format (construction via stream).
    template <typename format_type>
    void init_by_format()
    {
        static_assert(list_traits::contains<format_type, valid_formats>,
                      "You selected a format that is not in the valid_formats of this file.");

        format = detail::sam_file_input_format_exposer<format_type>{};
        secondary_stream = detail::make_secondary_istream(*primary_stream);
    }

    //!\brief The file header object.
    std::unique_ptr<header_type> header_ptr{new header_type{}};

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
                                         std::function<void(std::basic_istream<stream_char_type> *)>>;
    //!\brief Stream deleter that does nothing (no ownership assumed).
    static void stream_deleter_noop(std::basic_istream<stream_char_type> *)
    {}
    //!\brief Stream deleter with default behaviour (ownership assumed).
    static void stream_deleter_default(std::basic_istream<stream_char_type> * ptr)
    {
        delete ptr;
    }

    //!\brief The primary stream is the user provided stream or the file stream if constructed from filename.
    stream_ptr_t primary_stream{nullptr, stream_deleter_noop};
    //!\brief The secondary stream is a compression layer on the primary or just points to the primary (no compression).
    stream_ptr_t secondary_stream{nullptr, stream_deleter_noop};

    //!\brief Tracks whether the very first record is buffered when calling begin().
    bool first_record_was_read{false};
    //!\brief File is one position behind the last record.
    bool at_end{false};

    //!\brief Type of the format, a std::variant over the `valid_formats`.
    using format_type = typename detail::variant_from_tags<valid_formats, detail::sam_file_input_format_exposer>::type;

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
     * The reference information given by the IDs (`ref_ids`) and sequences (`ref_sequences`) will be used to keep the
     * record entry `seqan3::sam_file_input::record_type::reference_id()` consistent with the order imposed by
     * `ref_ids`. This way, you can use the value of `seqan3::sam_file_input::record_type::reference_id()` to access
     * the lists `ref_ids` and `ref_sequences` to retrieve the correct information for the current record.
     *
     * Note that the reference ids (names) must correspond to the exact spelling
     * in the SAM/BAM file. Otherwise, an exception will be thrown when reading.
     */
    template <std::ranges::forward_range ref_sequences_t>
    void set_references(typename traits_type::ref_ids & ref_ids, ref_sequences_t && ref_sequences)
    {
        assert(std::ranges::distance(ref_ids) == std::ranges::distance(ref_sequences));

        header_ptr = std::unique_ptr<header_type>{std::make_unique<header_type>(ref_ids)};
        reference_sequences_ptr = &ref_sequences;

        // initialise reference map and ref_dict if ref_ids are non-empty
        for (int32_t idx = 0; idx < std::ranges::distance(ref_ids); ++idx)
        {
            header_ptr->ref_id_info.emplace_back(std::ranges::distance(ref_sequences[idx]), "");

            if constexpr (std::ranges::contiguous_range<std::ranges::range_reference_t<typename traits_type::ref_ids>>
                          && std::ranges::sized_range<std::ranges::range_reference_t<typename traits_type::ref_ids>>
                          && std::ranges::borrowed_range<std::ranges::range_reference_t<typename traits_type::ref_ids>>)
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
        detail::get_or_ignore<field::header_ptr>(record_buffer) = header_ptr.get();

        // at end if we could not read further
        if (std::istreambuf_iterator<stream_char_type>{*secondary_stream}
            == std::istreambuf_iterator<stream_char_type>{})
        {
            at_end = true;
            return;
        }

        auto call_read_func = [this](auto & ref_seq_info)
        {
            std::visit(
                [&](auto & f)
                {
                    f.read_alignment_record(*secondary_stream,
                                            options,
                                            ref_seq_info,
                                            *header_ptr,
                                            position_buffer,
                                            detail::get_or_ignore<field::seq>(record_buffer),
                                            detail::get_or_ignore<field::qual>(record_buffer),
                                            detail::get_or_ignore<field::id>(record_buffer),
                                            detail::get_or_ignore<field::ref_seq>(record_buffer),
                                            detail::get_or_ignore<field::ref_id>(record_buffer),
                                            detail::get_or_ignore<field::ref_offset>(record_buffer),
                                            detail::get_or_ignore<field::cigar>(record_buffer),
                                            detail::get_or_ignore<field::flag>(record_buffer),
                                            detail::get_or_ignore<field::mapq>(record_buffer),
                                            detail::get_or_ignore<field::mate>(record_buffer),
                                            detail::get_or_ignore<field::tags>(record_buffer),
                                            detail::get_or_ignore<field::evalue>(record_buffer),
                                            detail::get_or_ignore<field::bit_score>(record_buffer));
                },
                format);
        };

        assert(!format.valueless_by_exception());

        if constexpr (!std::same_as<typename traits_type::ref_sequences, ref_info_not_given>)
            call_read_func(*reference_sequences_ptr);
        else
            call_read_func(std::ignore);
    }

    //!\brief Befriend iterator so it can access the buffers.
    friend iterator;
};

/*!\name Type deduction guides
 * \relates seqan3::sam_file_input
 * \{
 */
//!\brief Deduce selected fields, file_format, and default the rest.
template <input_stream stream_type, sam_file_input_format file_format, detail::fields_specialisation selected_field_ids>
sam_file_input(stream_type && stream, file_format const &, selected_field_ids const &)
    -> sam_file_input<typename sam_file_input<>::traits_type, // actually use the default
                      selected_field_ids,
                      type_list<file_format>>;

//!\brief Deduce selected fields, file_format, and default the rest.
template <input_stream stream_type, sam_file_input_format file_format, detail::fields_specialisation selected_field_ids>
sam_file_input(stream_type & stream, file_format const &, selected_field_ids const &)
    -> sam_file_input<typename sam_file_input<>::traits_type, // actually use the default
                      selected_field_ids,
                      type_list<file_format>>;

//!\brief Deduce file_format, and default the rest.
template <input_stream stream_type, sam_file_input_format file_format>
sam_file_input(stream_type && stream, file_format const &)
    -> sam_file_input<typename sam_file_input<>::traits_type,        // actually use the default
                      typename sam_file_input<>::selected_field_ids, // actually use the default
                      type_list<file_format>>;

//!\brief Deduce file_format, and default the rest.
template <input_stream stream_type, sam_file_input_format file_format>
sam_file_input(stream_type & stream, file_format const &)
    -> sam_file_input<typename sam_file_input<>::traits_type,        // actually use the default
                      typename sam_file_input<>::selected_field_ids, // actually use the default
                      type_list<file_format>>;

//!\brief Deduce selected fields, ref_sequences_t and ref_ids_t, default the rest.
template <std::ranges::forward_range ref_ids_t,
          std::ranges::forward_range ref_sequences_t,
          detail::fields_specialisation selected_field_ids>
sam_file_input(std::filesystem::path path, ref_ids_t &, ref_sequences_t &, selected_field_ids const &)
    -> sam_file_input<
        sam_file_input_default_traits<std::remove_reference_t<ref_sequences_t>, std::remove_reference_t<ref_ids_t>>,
        selected_field_ids,
        typename sam_file_input<>::valid_formats>; // actually use the default

//!\brief Deduce ref_sequences_t and ref_ids_t, default the rest.
template <std::ranges::forward_range ref_ids_t, std::ranges::forward_range ref_sequences_t>
sam_file_input(std::filesystem::path path, ref_ids_t &, ref_sequences_t &)
    -> sam_file_input<
        sam_file_input_default_traits<std::remove_reference_t<ref_sequences_t>, std::remove_reference_t<ref_ids_t>>,
        typename sam_file_input<>::selected_field_ids, // actually use the default
        typename sam_file_input<>::valid_formats>;     // actually use the default

//!\brief Deduce selected fields, ref_sequences_t and ref_ids_t, and file format.
template <input_stream stream_type,
          std::ranges::forward_range ref_ids_t,
          std::ranges::forward_range ref_sequences_t,
          sam_file_input_format file_format,
          detail::fields_specialisation selected_field_ids>
sam_file_input(stream_type && stream, ref_ids_t &, ref_sequences_t &, file_format const &, selected_field_ids const &)
    -> sam_file_input<
        sam_file_input_default_traits<std::remove_reference_t<ref_sequences_t>, std::remove_reference_t<ref_ids_t>>,
        selected_field_ids,
        type_list<file_format>>;

//!\brief Deduce selected fields, ref_sequences_t and ref_ids_t, and file format.
template <input_stream stream_type,
          std::ranges::forward_range ref_ids_t,
          std::ranges::forward_range ref_sequences_t,
          sam_file_input_format file_format,
          detail::fields_specialisation selected_field_ids>
sam_file_input(stream_type & stream, ref_ids_t &, ref_sequences_t &, file_format const &, selected_field_ids const &)
    -> sam_file_input<
        sam_file_input_default_traits<std::remove_reference_t<ref_sequences_t>, std::remove_reference_t<ref_ids_t>>,
        selected_field_ids,
        type_list<file_format>>;

//!\brief Deduce ref_sequences_t and ref_ids_t, and file format.
template <input_stream stream_type,
          std::ranges::forward_range ref_ids_t,
          std::ranges::forward_range ref_sequences_t,
          sam_file_input_format file_format>
sam_file_input(stream_type && stream, ref_ids_t &, ref_sequences_t &, file_format const &)
    -> sam_file_input<
        sam_file_input_default_traits<std::remove_reference_t<ref_sequences_t>, std::remove_reference_t<ref_ids_t>>,
        typename sam_file_input<>::selected_field_ids, // actually use the default
        type_list<file_format>>;

//!\brief Deduce selected fields, ref_sequences_t and ref_ids_t, and file format.
template <input_stream stream_type,
          std::ranges::forward_range ref_ids_t,
          std::ranges::forward_range ref_sequences_t,
          sam_file_input_format file_format>
sam_file_input(stream_type & stream, ref_ids_t &, ref_sequences_t &, file_format const &)
    -> sam_file_input<
        sam_file_input_default_traits<std::remove_reference_t<ref_sequences_t>, std::remove_reference_t<ref_ids_t>>,
        typename sam_file_input<>::selected_field_ids, // actually use the default
        type_list<file_format>>;
//!\}

} // namespace seqan3

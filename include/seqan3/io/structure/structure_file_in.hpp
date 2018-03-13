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
 * \brief Provides seqan3::structure_file_in and corresponding traits classes.
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 */

#pragma once

#include <cassert>
#include <fstream>
#include <limits>
#include <optional>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include <range/v3/algorithm/equal.hpp>

#include <seqan3/alphabet/adaptation/char.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
//#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/alphabet/nucleotide/rna5.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/structure/all.hpp>
#include <seqan3/core/metafunction/basic.hpp>
#include <seqan3/io/concept.hpp>
#include <seqan3/io/exception.hpp>
#include <seqan3/io/filesystem.hpp>
#include <seqan3/io/record.hpp>
#include <seqan3/io/detail/in_file_iterator.hpp>
#include <seqan3/io/detail/record.hpp>
#include <seqan3/io/structure/structure_file_in_format_concept.hpp>
#include <seqan3/io/structure/structure_file_format_dot_bracket.hpp>
#include <seqan3/range/container/concatenated_sequences.hpp>

namespace seqan3
{

//TODO(joergi-w) documentation
// ----------------------------------------------------------------------------
// sequence_file_in_traits_concept
// ----------------------------------------------------------------------------

/*!\interface seqan3::sequence_file_in_traits_concept <>
 * \brief The requirements a traits_type for seqan3::sequence_file_in must meet.
 * \ingroup sequence
 */
/*!\name Requirements for seqan3::sequence_file_in_traits_concept
 * \brief You can expect these **member types** of all types that satisfy seqan3::sequence_file_in_traits_concept.
 * \memberof seqan3::sequence_file_in_traits_concept
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
 * \memberof seqan3::sequence_file_in_traits_concept
 * \brief Alphabet of the characters for the seqan3::field::SEQ; must satisfy seqan3::alphabet_concept.
 */
/*!\typedef using sequence_legal_alphabet
 * \memberof seqan3::sequence_file_in_traits_concept
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
 * \memberof seqan3::sequence_file_in_traits_concept
 * \brief Type template of the seqan3::field::SEQ, a container template over `sequence_alphabet`;
 * must satisfy seqan3::sequence_concept.
 */
/*!\typedef using sequence_container_container
 * \memberof seqan3::sequence_file_in_traits_concept
 * \brief Type template of a column of seqan3::field::SEQ, a container template that can hold multiple
 * `sequence_container`; must satisfy seqan3::sequence_concept.
 */
/*!\typedef using id_alphabet
 * \memberof seqan3::sequence_file_in_traits_concept
 * \brief Alphabet of the characters for the seqan3::field::ID; must satisfy seqan3::alphabet_concept.
 */
/*!\typedef using id_container
 * \memberof seqan3::sequence_file_in_traits_concept
 * \brief Type template of the seqan3::field::ID, a container template over `id_alphabet`;
 * must satisfy seqan3::sequence_concept.
 */
/*!\typedef using id_container_container
 * \memberof seqan3::sequence_file_in_traits_concept
 * \brief Type template of a column of seqan3::field::ID, a container template that can hold multiple
 * `id_container`; must satisfy seqan3::sequence_concept.
 */
/*!\typedef using quality_alphabet
 * \memberof seqan3::sequence_file_in_traits_concept
 * \brief Alphabet of the characters for the seqan3::field::QUAL; must satisfy seqan3::quality_concept.
 */
/*!\typedef using quality_container
 * \memberof seqan3::sequence_file_in_traits_concept
 * \brief Type template of the seqan3::field::QUAL, a container template over `quality_alphabet`;
 * must satisfy seqan3::sequence_concept.
 */
/*!\typedef using quality_container_container
 * \memberof seqan3::sequence_file_in_traits_concept
 * \brief Type template of a column of seqan3::field::QUAL, a container template that can hold multiple
 * `quality_container`; must satisfy seqan3::sequence_concept.
 */
//!\}
//!\cond
template<typename t> concept bool structure_file_in_traits_concept = requires(t v)
{
    // sequence
    requires alphabet_concept<typename t::seq_alphabet>;
    requires alphabet_concept<typename t::seq_legal_alphabet>;
    requires explicitly_convertible_to_concept<typename t::seq_legal_alphabet, typename t::seq_alphabet>;
    requires sequence_container_concept<typename t::template seq_container<typename t::seq_alphabet>>;
    requires sequence_container_concept
        <typename t::template seq_container_container
            <typename t::template seq_container
                <typename t::seq_alphabet>>>;

    // id
    requires alphabet_concept<typename t::id_alphabet>;
    requires sequence_container_concept<typename t::template id_container<typename t::id_alphabet>>;
    requires sequence_container_concept
        <typename t::template id_container_container
            <typename t::template id_container
                <typename t::id_alphabet>>>;

    // bpp
    requires std::is_floating_point_v<typename t::bpp_prec>;
    requires std::numeric_limits<typename t::bpp_partner>::is_integer;
    requires sequence_container_concept
        <typename t::template bpp_container
            <typename t::template bpp_queue
                 <typename t::template bpp_item
                      <typename t::bpp_prec, typename t::bpp_partner>>>>;
    requires sequence_container_concept
        <typename t::template bpp_container_container
            <typename t::template bpp_container
                <typename t::template bpp_queue
                    <typename t::template bpp_item
                        <typename t::bpp_prec, typename t::bpp_partner>>>>>;

    // structure
    requires rna_structure_concept<typename t::structure_alphabet>;
    requires rna_structure_concept<typename t::structure_legal_alphabet>;
    requires explicitly_convertible_to_concept<typename t::structure_legal_alphabet, typename t::structure_alphabet>;
    requires sequence_container_concept<typename t::template structure_container<typename t::structure_alphabet>>;
    requires sequence_container_concept
        <typename t::template structure_container_container
            <typename t::template structure_container
                <typename t::structure_alphabet>>>;

    // structured sequence
    requires sequence_container_concept
        <typename t::template structured_seq_container
            <typename t::template structured_seq_alphabet
                <typename t::seq_alphabet, typename t::structure_alphabet>>>;
    requires sequence_container_concept
        <typename t::template structured_seq_container_container
            <typename t::template structured_seq_container
                <typename t::template structured_seq_alphabet
                    <typename t::seq_alphabet, typename t::structure_alphabet>>>>;

    // energy: std::optional of floating point number
    requires std::is_floating_point_v<typename t::energy_type::value_type>;
    requires sequence_container_concept<typename t::template energy_container<typename t::energy_type>>;

    // reactivity [error]
    requires std::is_floating_point_v<typename t::react_type>;
    requires sequence_container_concept<typename t::template react_container<typename t::react_type>>;
    requires sequence_container_concept
        <typename t::template react_container_container
            <typename t::template react_container
                <typename t::react_type>>>;

    // comment
    requires alphabet_concept<typename t::comment_alphabet>;
    requires sequence_container_concept<typename t::template comment_container<typename t::comment_alphabet>>;
    requires sequence_container_concept
        <typename t::template comment_container_container
            <typename t::template comment_container
                <typename t::comment_alphabet>>>;

    // offset
    requires std::numeric_limits<typename t::offset_type>::is_integer;
    requires sequence_container_concept<typename t::template offset_container<typename t::offset_type>>;
};
//!\endcond

// ----------------------------------------------------------------------------
// structure_file_in_default_traits
// ----------------------------------------------------------------------------

/*!\brief The default traits for seqan3::structure_file_in
 * \implements structure_file_in_traits_concept
 * \ingroup structure
 *
 * \details
 *
 * If you wish to change a single or a few types from the default, just inherit from this class and
 * "overwrite" the respective type definitions.
 *
 * This example will make the file read into a smaller alphabet and a compressed container:
 *
 * ```cpp
 * struct my_traits : structure_file_in_default_traits_dna
 * {
 *     using sequence_alphabet = dna4;                        // instead of dna5
 *
 *     template <typename alph>
 *     using sequence_container = bitcompressed_vector<alph>; // must be defined as a template!
 * };
 *
 * structure_file_in<my_traits> fin{"/tmp/my.fasta"};
 *
 * //...
 * ```
 */
struct structure_file_in_default_traits_rna
{
    /*!\name Member types
     * \brief Definitions to satisfy seqan3::structure_file_in_traits_concept.
     * \{
     */

    // sequence
    using seq_alphabet                       = rna5;
    using seq_legal_alphabet                 = rna15;
    template<typename _seq_alphabet>
    using seq_container                      = std::vector<_seq_alphabet>;
    template<typename _seq_container>
    using seq_container_container            = concatenated_sequences<_seq_container>;

    // id
    using id_alphabet                        = char;
    template<typename _id_alphabet>
    using id_container                       = std::basic_string<_id_alphabet>;
    template<typename _id_container>
    using id_container_container             = concatenated_sequences<_id_container>;

    // base pair probability structure
    using bpp_prec                           = double;
    using bpp_partner                        = size_t;
    template<typename _bpp_prec, typename _bpp_partner>
    using bpp_item                           = std::pair<_bpp_prec, _bpp_partner>;
    template<typename _bpp_item>
    using bpp_queue                          = std::set<_bpp_item>;
    template<typename _bpp_queue>
    using bpp_container                      = std::vector<_bpp_queue>;
    template<typename _bpp_container>
    using bpp_container_container            = std::vector<_bpp_container>;

    // fixed structure
    using structure_alphabet                 = wuss51;
    using structure_legal_alphabet           = wuss51;
    template<typename _structure_alphabet>
    using structure_container                = std::vector<_structure_alphabet>;
    template<typename _structure_container>
    using structure_container_container      = concatenated_sequences<_structure_container>;

    // combined sequence and structure
    template<typename _sequence_alphabet, typename _structure_alphabet>
    using structured_seq_alphabet            = structured_rna<_sequence_alphabet, _structure_alphabet>;
    template<typename _structured_seq_alphabet>
    using structured_seq_container           = std::vector<_structured_seq_alphabet>;
    template<typename _structured_seq_container>
    using structured_seq_container_container = concatenated_sequences<_structured_seq_container>;

    // energy
    using energy_type                        = std::optional<double>;
    template<typename _energy_type>
    using energy_container                   = std::vector<_energy_type>;

    // reactivity [error]
    using react_type                         = double;
    template<typename _react_type>
    using react_container                    = std::vector<_react_type>;
    template<typename _react_container>
    using react_container_container          = std::vector<_react_container>;

    // comment
    using comment_alphabet                   = char;
    template<typename _comment_alphabet>
    using comment_container                  = std::basic_string<_comment_alphabet>;
    template<typename _comment_container>
    using comment_container_container        = concatenated_sequences<_comment_container>;

    // offset
    using offset_type                        = size_t;
    template<typename _offset_type>
    using offset_container                   = std::vector<_offset_type>;
    //!\}
};

struct structure_file_in_default_traits_aa : structure_file_in_default_traits_rna
{
    /*!\name Member types
     * \brief Definitions to satisfy seqan3::structure_file_in_traits_concept.
     * \{
     */
    using seq_alphabet        = aa27;
    using seq_legal_alphabet  = aa27;
    using structure_alphabet       = dssp9;
    using structure_legal_alphabet = dssp9;
    template<typename _seq_alphabet, typename _structure_alphabet>
    using structured_seq_alphabet        = structured_aa<_seq_alphabet, _structure_alphabet>;
    //!\}
};

// ----------------------------------------------------------------------------
// structure_file_in
// ----------------------------------------------------------------------------

/*!\brief A class for reading sequence files, e.g. FASTA, FASTQ ...
 * \ingroup io
 * \tparam traits_type An auxiliary type that defines certain member type and constants, must satisfy
 * seqan3::sequence_file_traits_concept.
 * \tparam stream_type The type of the stream, must satisfy seqan3::istream_concept.
 * \tparam selected_fields A \ref fields type with the list and order of desired record entries.
 *
 * \details
 *
 * ### Introduction
 *
 * Sequence files are the most generic and common biological files. Well-known formats include
 * FastA and FastQ, but some may also be interested in treating SAM or BAM files as sequence
 * files, discarding the alignment.
 *
 * The Sequence file abstraction provides two fields: seqan3::field::SEQ and seqan3::field::ID.
 *
 * By default the SEQ field is retrieved as a vector over seqan3::qualified <seqan3::dna5>, i.e.
 * the sequence combines qualities and actual sequence in one. You can later drop the qualities if
 * you are not interested in them, or specify a custom traits type to change the underlying
 * alphabet of the sequence field so they are never returned.
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
 * structure_file_in sf{"/tmp/my.db"}; // Dot bracket with RNA sequences assumed, regular std::ifstream taken as stream
 * ```
 * Reading from a std::istringstream:
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
 * sequence_file_in fin{std::move(iss), sequence_file_format_fasta{}};
 * //              ^ no need to specify the template arguments
 * ```
 *
 * In some cases, you do need to specify the arguments, e.g. if you want to read amino acids:
 *
 * ```cpp
 * sequence_file_in<seqeuence_file_default_traits_aa> fin{"/tmp/my.fasta"};
 * ```
 *
 * You can define your own traits type to further customise the types used by and return by this class, see
 * seqan3::seqeuence_file_default_traits_dna for more details.
 *
 * ### Reading record-wise
 *
 * You can iterate over this file record-wise:
 *
 * ```cpp
 * sequence_file_in fin{"/tmp/my.fasta"};
 *
 * for (auto record : fin)
 * {
 *     std::cout << "ID:  " << std::get<1>(record) << '\n';
 *     std::cout << "SEQ: " << (std::get<0>(record) | view::to_char) << '\n'; // sequence is converted to char string on-the-fly
 * }
 * ```
 *
 * In the above example, record has the type \ref record_type which is a tuple, that's why we can access it via
 * std::get.
 *
 * However we can also directly use [structured bindings](http://en.cppreference.com/w/cpp/language/structured_binding)
 * to decompose the tuple into its elements:
 *
 * ```cpp
 * sequence_file_in fin{"/tmp/my.fasta"};
 *
 * for (auto [ seq, id ] : fin)
 * {
 *     std::cout << "ID:  " << id << '\n';
 *     std::cout << "SEQ: " << (seq | view::to_char) << '\n'; // sequence string is converted to char string on-the-fly
 * }
 * ```
 * In this case you immediately get the two elements of the tuple: `seq` of \ref sequence_type and `id` of
 * \ref id_type.
 *
 * ### Reading column-wise
 *
 * TODO example
 *
 *
 * ### Formats
 *
 * TODO give overview of formats
 *
 *
 *
 */
template<structure_file_in_traits_concept traits_type_ = structure_file_in_default_traits_rna,
         detail::fields_concept selected_field_ids_ = fields<field::SEQ, field::ID, field::BPP>,
         detail::type_list_of_structure_file_in_formats_concept valid_formats_
             = type_list<structure_file_format_dot_bracket>,
         istream_concept<char> stream_type_ = std::ifstream>
class structure_file_in
{
public:
    /*!\name Template arguments
     * \brief Exposed as member types for public access.
     * \{
     */
    //!\brief A traits type that defines aliases and template for storage of the fields.
    using traits_type        = traits_type_;
    //!\brief A seqan3::fields list with the fields selected for the record.
    using selected_field_ids = selected_field_ids_;
    //!\brief A seqan3::type_list with the possible formats.
    using valid_formats      = valid_formats_;
    //!\brief The type of the underlying stream.
    using stream_type        = stream_type_;
    //!\}

    /*!\brief The subset of seqan3::field IDs that are valid for this file; order corresponds to the types in
     * \ref field_types.
      */
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

    static_assert([]() constexpr
                  {
                      for (field f : selected_field_ids::as_array)
                          if (!field_ids::contains(f))
                              return false;
                      return true;
                  }(),
                  "You selected a field that is not valid for structure files, please refer to the documentation "
                  "of structure_file_in::field_ids for the accepted values.");

    static_assert([]() constexpr
                  {
                      return !(selected_field_ids::contains(field::STRUCTURED_SEQ) &&
                               (selected_field_ids::contains(field::SEQ) ||
                               (selected_field_ids::contains(field::STRUCTURE))));
                  }(), "You may not select field::STRUCTURED_SEQ and either of field::SEQ and field::STRUCTURE "
                       "at the same time.");

    /*!\name Field types and record type
     * \brief These types are relevant for record/row-based reading; they may be manipulated via the \ref traits_type
     * to achieve different storage behaviour.
     * \{
     */
    //!\brief The type of the sequence field (default std::vector of seqan3::rna5).
    using seq_type            = typename traits_type::template seq_container<typename traits_type::seq_alphabet>;
    //!\brief The type of the ID field (default std::string).
    using id_type             = typename traits_type::template id_container<typename traits_type::id_alphabet>;
    //!\brief The type of the base pair probabilies (default std::vector of std::set<std::pair<double, size_t>>).
    using bpp_type            = typename traits_type::template bpp_container
        <typename traits_type::template bpp_queue
            <typename traits_type::template bpp_item
                <typename traits_type::bpp_partner, typename traits_type::bpp_prec>>>;
    //!\brief The type of the structure field (default std::vector of seqan3::wuss51).
    using structure_type      = typename traits_type::template structure_container
        <typename traits_type::structure_alphabet>;
    //!\brief The type of the sequence-structure field (default std::vector of structured_rna<rna5,wuss51>).
    using structured_seq_type = typename traits_type::template structured_seq_container
        <typename traits_type::template structured_seq_alphabet
            <typename traits_type::seq_alphabet, typename traits_type::structure_alphabet>>;
    //!\brief The type of the energy field (default double).
    using energy_type         = typename traits_type::energy_type;
    //!\brief The type of the energy field (default std::string).
    using react_type          = typename traits_type::template react_container<typename traits_type::react_type>;
    //!\brief The type of the energy field (default double).
    using comment_type        = typename traits_type::template comment_container
        <typename traits_type::comment_alphabet>;
    //!\brief The type of the offset field (default size_t).
    using offset_type         = typename traits_type::offset_type;

    /* ALIGN */

    //!\brief The previously defined types in a type list (meta::list).
    using field_types    = type_list<seq_type, id_type, bpp_type, structure_type, structured_seq_type, energy_type,
                                     react_type, react_type, comment_type, offset_type>;

    //!\brief The type of the record, specialisation of seqan3::record; acts as a tuple of the selected field types.
    using record_type    = record<detail::select_types_with_ids_t<field_types, field_ids, selected_field_ids>,
                                  selected_field_ids>;
    //!\}

    /*!\name Field column types and tuple type
     * \brief These types are relevant for field/column-wise reading; they may be manipulated via the \ref traits_type
     * to achieve different storage behaviour.
     * \{
     */
     //!\brief Column type of field::SEQ (seqan3::concatenated_sequences<seq_type> by default).
    using seq_column_type            = typename traits_type::template seq_container_container<seq_type>;
    //!\brief Column type of field::ID (seqan3::concatenated_sequences<id_type> by default).
    using id_column_type             = typename traits_type::template id_container_container<id_type>;
    //!\brief Column type of field::BPP (std::vector<bpp_type> by default).
    using bpp_column_type            = typename traits_type::template bpp_container_container<bpp_type>;
    //!\brief Column type of field::STRUCTURE (seqan3::concatenated_sequences<structure_type> by default).
    using structure_column_type      = typename traits_type::template structure_container_container<structure_type>;
    //!\brief Column type of field::STRUCTURED_SEQ (seqan3::concatenated_sequences<structured_seq_type> by default).
    using structured_seq_column_type = typename traits_type::template structured_seq_container_container
        <structured_seq_type>;
    //!\brief Column type of field::ENERGY (std::vector<energy_type> by default).
    using energy_column_type         = typename traits_type::template energy_container<energy_type>;
    //!\brief Column type of field::REACT and field::REACT_ERR (std::vector<react_type> by default).
    using react_column_type          = typename traits_type::template react_container_container<react_type>;
    //!\brief Column type of field::COMMENT (seqan3::concatenated_sequences<comment_type> by default).
    using comment_column_type        = typename traits_type::template comment_container_container<comment_type>;
    //!\brief Column type of field::OFFSET (std::vector<offset_type> by default).
    using offset_column_type         = typename traits_type::template offset_container<offset_type>;

    //!\brief The previously defined types aggregated in a seqan3::type_list.
    using field_column_types         = type_list<seq_column_type,
                                                 id_column_type,
                                                 bpp_column_type,
                                                 structure_column_type,
                                                 structured_seq_column_type,
                                                 energy_column_type,
                                                 react_column_type,
                                                 react_column_type,
                                                 comment_column_type,
                                                 offset_column_type>;
    //!\brief The type emulated by the file when read column-wise.
    using file_as_tuple_type     = record<detail::select_types_with_ids_t<field_column_types, field_ids,
                                                                          selected_field_ids>, selected_field_ids>;
    //!\}

    /*!\name Range associated types
     * \brief The types necessary to facilitate the behaviour of an input range (used in record-wise reading).
     * \{
     */
    //!\brief The value_type is the \ref record_type.
    using value_type      = record_type;
    //!\brief The reference type.
    using reference       = record_type &;
    //!\brief The const_reference type is void, because files are not const-iterable.
    using const_reference = void;
    //!\brief An unsigned integer type, usually std::size_t.
    using size_type       = size_t;
    //!\brief A signed integer type, usually std::ptrdiff_t.
    using difference_type = std::make_signed_t<size_t>;
    //!\brief The iterator type of this view (an input iterator).
    using iterator        = detail::in_file_iterator<structure_file_in>;
    //!\brief The const iterator type is void, because files are not const-iterable.
    using const_iterator  = void;
    //!\brief The type returned by end().
    using sentinel        = ranges::default_sentinel;
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    //!\brief Default constructor is explicitly deleted, you need to give a stream or file name.
    structure_file_in() = delete;
    //!\brief Copy construction is explicitly deleted, because you cannot have multiple access to the same file.
    structure_file_in(structure_file_in const &) = delete;
    //!\brief Copy assignment is explicitly deleted, because you cannot have multiple access to the same file.
    structure_file_in & operator=(structure_file_in const &) = delete;
    //!\brief Move construction is defaulted.
    structure_file_in(structure_file_in &&) = default;
    //!\brief Move assignment is defaulted.
    structure_file_in & operator=(structure_file_in &&) = default;
    //!\brief Destructor is defaulted.
    ~structure_file_in() = default;

    /*!\brief Construct from filename.
     * \param[in] _file_name    Path to the file you wish to open.
     * \param[in] fields_tag    A seqan3::fields tag. [optional]
     *
     * \details
     *
     * In addition to the file name, you may specify a custom seqan3::fields type which may be easier than
     * defining all the template parameters.
     */
    structure_file_in(filesystem::path const & _file_name,
                      selected_field_ids const & SEQAN3_DOXYGEN_ONLY(fields_tag) = selected_field_ids{})
    {
        // open stream
        stream.open(_file_name, std::ios_base::in | std::ios::binary);
        if (!stream.is_open())
            throw file_open_error{"Could not open file for reading."};

        // initialise format handler
        bool format_found = false;
        std::string extension = _file_name.extension().string();
        if (extension.size() > 1)
        {
            extension = extension.substr(1); // drop leading "."
            meta::for_each(valid_formats{}, [&] (auto && fmt)
            {
                using fmt_type = remove_cvref_t<decltype(fmt)>;

                for (auto const & ext : fmt_type::file_extensions)
                {
                    if (ranges::equal(ext, extension))
                    {
                        format = fmt_type{};
                        format_found = true;
                        return;
                    }
                }
            });
        }
        if (!format_found)
            throw unhandled_extension_error("No valid format found for this extension.");

        // buffer first record
        read_next_record();
    }

    /*!\brief Construct from an existing stream and with specified format.
     * \tparam file_format   The format of the file in the stream, must satisfy seqan3::structure_file_in_format_concept.
     * \param[in] _stream    The stream to operate on (this must be std::move'd in!).
     * \param[in] format_tag The file format tag.
     * \param[in] fields_tag A seqan3::fields tag. [optional]
     */
    template<structure_file_in_format_concept file_format>
    structure_file_in(stream_type && _stream,
                      file_format const & SEQAN3_DOXYGEN_ONLY(format_tag),
                      selected_field_ids const & SEQAN3_DOXYGEN_ONLY(fields_tag) = selected_field_ids{}) :
        stream{std::move(_stream)}, format{file_format{}}
    {
        static_assert(meta::in<valid_formats, file_format>::value,
                      "You selected a format that is not in the valid_formats of this file.");

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
     * structure_file_in fin{"/tmp/my.dbn"};
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
     * sequence_file_in fin{"/tmp/my.dbn"};
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
    template<field f>
    friend auto & get(structure_file_in & file)
    {
        static_assert(structure_file_in::selected_field_ids::contains(f),
                      "You requested a field via get that was not selected for the file.");

        file.read_columns();

        return seqan3::get<f>(file.columns_buffer);
    }

    //!\copydoc get
    template<field f>
    friend auto && get(structure_file_in && file)
    {
        return std::move(get<f>(file));
    }

    //!\copydoc get
    template<size_t i>
    friend auto & get(structure_file_in & file)
    {
        static_assert(i < structure_file_in::selected_field_ids::as_array.size(),
                      "You requested a field number larger than the number of selected fields for the file.");
        file.read_columns();

        return std::get<i>(file.columns_buffer);
    }

    //!\copydoc get
    template<size_t i>
    friend auto && get(structure_file_in && file)
    {
        return std::move(get<i>(file));
    }

    //!\copydoc get
    template<typename t>
    friend auto & get(structure_file_in & file)
    {
        file.read_columns();

        return std::get<t>(file.columns_buffer);
    }

    //!\copydoc get
    template<typename t>
    friend auto && get(structure_file_in && file)
    {
        return std::move(get<t>(file));
    }
    //!\}

    //!\brief The options are public and its members can be set directly.
    structure_file_in_options<typename traits_type::seq_legal_alphabet,
                              typename traits_type::structure_legal_alphabet> options;

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

    //!\brief Path of the file that the stream operates on.
    std::string file_name;

    //!\brief The stream we are reading from.
    stream_type stream;

    //!\brief File is at position 1 behind the last record.
    bool at_end{false};

    //!\brief Type of the format, an std::variant over `valid_formats`.
    using format_type = detail::transfer_template_args_onto_t<valid_formats, std::variant>;
    //!\brief The actual std::variant holding a pointer to the detected/selected format.
    format_type format;

    //!\brief Tell the format to move to the next record and update the buffer.
    void read_next_record()
    {
        if (at_end)
            return;

        // clear the record
        record_buffer.clear();

        // at end if we could not read further
        if (stream.eof())
        {
            at_end = true;
            return;
        }

        assert(!format.valueless_by_exception());
        std::visit([&] (structure_file_in_format_concept & f)
        {
            // read new record
            f.read(stream,
                   options,
                   detail::get_or_ignore<field::SEQ>(record_buffer),
                   detail::get_or_ignore<field::ID>(record_buffer),
                   detail::get_or_ignore<field::BPP>(record_buffer),
                   detail::get_or_ignore<field::STRUCTURE>(record_buffer),
                   detail::get_or_ignore<field::STRUCTURED_SEQ>(record_buffer),
                   detail::get_or_ignore<field::ENERGY>(record_buffer),
                   detail::get_or_ignore<field::REACT>(record_buffer),
                   detail::get_or_ignore<field::REACT_ERR>(record_buffer),
                   detail::get_or_ignore<field::COMMENT>(record_buffer),
                   detail::get_or_ignore<field::OFFSET>(record_buffer));
        }, format);
    }

    //!\brief Read the entire file into the internal column buffers.
    void read_columns()
    {
//        std::cout << "read_columns" << std::endl;
        auto & seq_column_buffer            = detail::get_or_ignore<field::SEQ>(columns_buffer);
        auto & id_column_buffer             = detail::get_or_ignore<field::ID>(columns_buffer);
        auto & bpp_column_buffer            = detail::get_or_ignore<field::BPP>(columns_buffer);
        auto & structure_column_buffer      = detail::get_or_ignore<field::STRUCTURE>(columns_buffer);
        auto & structured_seq_column_buffer = detail::get_or_ignore<field::STRUCTURED_SEQ>(columns_buffer);
        auto & energy_column_buffer         = detail::get_or_ignore<field::ENERGY>(columns_buffer);
        auto & react_column_buffer          = detail::get_or_ignore<field::REACT>(columns_buffer);
        auto & react_err_column_buffer      = detail::get_or_ignore<field::REACT_ERR>(columns_buffer);
        auto & comment_column_buffer        = detail::get_or_ignore<field::COMMENT>(columns_buffer);
        auto & offset_column_buffer         = detail::get_or_ignore<field::OFFSET>(columns_buffer);

        if (true || !stream.eof())
        {
            // read the remaining records and split into column buffers
            for (auto & rec : *this)
            {
                if constexpr (selected_field_ids::contains(field::SEQ))
                    seq_column_buffer.push_back(std::move(seqan3::get<field::SEQ>(rec)));
                if constexpr (selected_field_ids::contains(field::ID))
                    id_column_buffer.push_back(std::move(seqan3::get<field::ID>(rec)));
                if constexpr (selected_field_ids::contains(field::BPP))
                    bpp_column_buffer.push_back(std::move(seqan3::get<field::BPP>(rec)));
                if constexpr (selected_field_ids::contains(field::STRUCTURE))
                    structure_column_buffer.push_back(std::move(seqan3::get<field::STRUCTURE>(rec)));
                if constexpr (selected_field_ids::contains(field::STRUCTURED_SEQ))
                    structured_seq_column_buffer.push_back(std::move(seqan3::get<field::STRUCTURED_SEQ>(rec)));
                if constexpr (selected_field_ids::contains(field::ENERGY))
                    energy_column_buffer.push_back(std::move(seqan3::get<field::ENERGY>(rec)));
                if constexpr (selected_field_ids::contains(field::REACT))
                    react_column_buffer.push_back(std::move(seqan3::get<field::REACT>(rec)));
                if constexpr (selected_field_ids::contains(field::REACT_ERR))
                    react_err_column_buffer.push_back(std::move(seqan3::get<field::REACT_ERR>(rec)));
                if constexpr (selected_field_ids::contains(field::COMMENT))
                    comment_column_buffer.push_back(std::move(seqan3::get<field::COMMENT>(rec)));
                if constexpr (selected_field_ids::contains(field::OFFSET))
                    offset_column_buffer.push_back(std::move(seqan3::get<field::OFFSET>(rec)));
            }
        }
    }

    //!\brief Befriend iterator so it can access the buffers.
    friend iterator;
};

/*!\name Type deduction guides
 * \relates seqan3::structure_file_in
 * \{
 */
template <istream_concept<char> stream_type,
          structure_file_in_format_concept file_format,
          detail::fields_concept selected_field_ids>
structure_file_in(stream_type && _stream, file_format const &, selected_field_ids const &)
    -> structure_file_in<typename structure_file_in<>::traits_type,       // actually use the default
                         selected_field_ids,
                         type_list<file_format>,
                         std::remove_reference_t<stream_type>>;
//!\}

} // namespace seqan3

// ------------------------------------------------------------------
// std-overloads for the tuple-like interface
// ------------------------------------------------------------------

namespace std
{
//!\brief std::tuple_size overload for column-like access. [metafunction specialisation for seqan3::structure_file_in]
template<seqan3::structure_file_in_traits_concept traits_type,
         seqan3::detail::fields_concept selected_field_ids,
         seqan3::detail::type_list_of_structure_file_in_formats_concept valid_formats,
         seqan3::istream_concept<char> stream_type>
struct tuple_size<seqan3::structure_file_in<traits_type, selected_field_ids, valid_formats, stream_type>>
{
    //!\brief The value equals the number of selected fields in the file.
    static constexpr size_t value = selected_field_ids::as_array.size();
};

//!\brief std::tuple_element overload for column-like access. [metafunction specialisation for seqan3::structure_file_in]
template<size_t elem_no,
         seqan3::structure_file_in_traits_concept traits_type,
         seqan3::detail::fields_concept selected_field_ids,
         seqan3::detail::type_list_of_structure_file_in_formats_concept valid_formats,
         seqan3::istream_concept<char> stream_type>
struct tuple_element<elem_no, seqan3::structure_file_in<traits_type, selected_field_ids, valid_formats, stream_type>>
    : tuple_element<elem_no,
                    typename seqan3::structure_file_in<traits_type,
                                                       selected_field_ids,
                                                       valid_formats,
                                                       stream_type>::file_as_tuple_type>
{};

} // namespace std

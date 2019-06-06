// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::structure_file_input and corresponding traits classes.
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
#include <seqan3/alphabet/aminoacid/all.hpp>
#include <seqan3/alphabet/nucleotide/rna15.hpp>
#include <seqan3/alphabet/nucleotide/rna5.hpp>
#include <seqan3/alphabet/structure/all.hpp>
#include <seqan3/core/type_traits/basic.hpp>
#include <seqan3/io/stream/concept.hpp>
#include <seqan3/io/exception.hpp>
#include <seqan3/std/filesystem>
#include <seqan3/io/record.hpp>
#include <seqan3/io/detail/in_file_iterator.hpp>
#include <seqan3/io/detail/misc_input.hpp>
#include <seqan3/io/detail/record.hpp>
#include <seqan3/io/structure_file/input_format_concept.hpp>
#include <seqan3/io/structure_file/input_options.hpp>
#include <seqan3/io/structure_file/format_vienna.hpp>
#include <seqan3/range/container/concatenated_sequences.hpp>

namespace seqan3
{
// ----------------------------------------------------------------------------
// StructureFileInputTraits
// ----------------------------------------------------------------------------

/*!\interface seqan3::StructureFileInputTraits <>
 * \brief The requirements a traits_type for seqan3::structure_file_input must meet.
 * \ingroup structure_file
 */
/*!\name Requirements for seqan3::StructureFileInputTraits
 * \brief You can expect these **member types** of all types that satisfy seqan3::StructureFileInputTraits.
 * \memberof seqan3::StructureFileInputTraits
 *
 * \details
 *
 * Note that the alphabet type of the seqan3::field::STRUCTURED_SEQ is a combined alphabet, i.e. either
 * seqan3::structured_rna<sequence_alphabet, structure_alphabet> or
 * seqan3::structured_aa<sequence_alphabet, structure_alphabet> and the container type templates for
 * the field are those of seqan3::field::SEQ and seqan3::field::STRUCTURE, respectively.
 *
 * \{
 */
/*!\typedef using seq_alphabet
 * \brief Alphabet of the characters for the seqan3::field::SEQ; must satisfy seqan3::Alphabet.
 */
/*!\typedef using seq_legal_alphabet
 * \brief Intermediate alphabet for seqan3::field::SEQ; must satisfy seqan3::Alphabet and be convertible to
 * `seq_alphabet`.
 *
 * \details
 *
 * This alphabet can be a superset of `seq_alphabet` to allow conversion of some characters
 * without producing an error, e.g. if this is set to seqan3::rna15 and `seq_alphabet` is set to seqan3::rna5,
 * 'M' will be an accepted character and automatically converted to 'N', while 'Z' will still be an illegal
 * character and produce an error.
 */
/*!\typedef using seq_container
 * \brief Type template of the seqan3::field::SEQ, a container template over `seq_alphabet`;
 * must satisfy seqan3::SequenceContainer.
 */
/*!\typedef using seq_container_container
 * \brief Type template of a column of seqan3::field::SEQ, a container template that can hold multiple
 * `seq_container`; must satisfy seqan3::SequenceContainer.
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
/*!\typedef using bpp_prob
 * \brief Data type for the base pair probabilities in seqan3::field::BPP; must satisfy std::is_floating_point.
 */
/*!\typedef using bpp_partner
 * \brief Data type for the partner index of an interaction in seqan3::field::BPP; must satisfy
 * std::numeric_limits::is_integer.
 */
/*!\typedef using bpp_queue
 * \brief A container template representing a set of interactions of type bpp_item,
 * which are (comparable) tuples of `bpp_prob` and `bpp_partner`;
 * must satisfy seqan3::Container and must provide an std::emplace(bpp_prob, bpp_partner) function.
 */
/*!\typedef using bpp_container
 * \brief Type template of the seqan3::field::BPP, a container template over a set (bpp_queue) of interactions;
 * must satisfy seqan3::SequenceContainer.
 */
/*!\typedef using bpp_container_container
 * \brief Type template of a column of seqan3::field::BPP, a container template that can hold multiple
 * `bpp_container`; must satisfy seqan3::SequenceContainer.
 */
/*!\typedef using structure_alphabet
 * \brief Alphabet of the characters for the seqan3::field::STRUCTURE; must satisfy seqan3::RnaStructureAlphabet.
 */
/*!\typedef using structure_container
 * \brief Type template of the seqan3::field::STRUCTURE, a container template over `structure_alphabet`;
 * must satisfy seqan3::SequenceContainer.
 */
/*!\typedef using structure_container_container
 * \brief Type template of a column of seqan3::field::STRUCTURE, a container template that can hold multiple
 * `structure_container`; must satisfy seqan3::SequenceContainer.
 */
/*!\typedef using energy_type
 * \brief Type template of the seqan3::field::ENERGY; must be std::optional of a type satisfying std::is_floating_point.
 * \details
 * If the file record contains an energy, the value can be obtained through dereference or std::optional::value.
 * Otherwise, operator bool or std::optional::has_value return false.
 */
/*!\typedef using energy_container
 * \brief Type template of a column of seqan3::field::ENERGY, a container template that can hold multiple `energy_type`;
 * must satisfy seqan3::SequenceContainer.
 */
/*!\typedef using react_type
 * \brief Data type for the reactivity and reactivity error in seqan3::field::REACT and seqan3::field::REACT_ERR,
 * respectively; must satisfy std::is_floating_point.
 */
/*!\typedef using react_container
 * \brief Type template of the seqan3::field::REACT and seqan3::field::REACT_ERR, a container template over
 * `react_type`; must satisfy seqan3::SequenceContainer.
 */
/*!\typedef using react_container_container
 * \brief Type template of a column of seqan3::field::REACT and seqan3::field::REACT_ERR, a container template that
 * can hold multiple `react_container`; must satisfy seqan3::SequenceContainer.
 */
/*!\typedef using comment_alphabet
 * \brief Alphabet of the characters for the seqan3::field::COMMENT; must satisfy seqan3::Alphabet.
 */
/*!\typedef using comment_container
 * \brief Type template of the seqan3::field::COMMENT, a container template over `comment_alphabet`;
 * must satisfy seqan3::SequenceContainer.
 */
/*!\typedef using comment_container_container
 * \brief Type template of a column of seqan3::field::COMMENT, a container template that can hold multiple
 * `comment_container`; must satisfy seqan3::SequenceContainer.
 */
/*!\typedef using offset_type
 * \brief Type template of the seqan3::field::OFFSET; must statisfy std::numeric_limits::is_integer.
 */
/*!\typedef using offset_container
 * \brief Type template of a column of seqan3::field::OFFSET, a container template that can hold multiple `offset_type`;
 * must satisfy seqan3::SequenceContainer.
 */
//!\}
//!\cond
template<typename t>
SEQAN3_CONCEPT StructureFileInputTraits = requires(t v)
{
    // TODO(joergi-w) The expensive concept checks are currently omitted. Check again when compiler has improved.
    // sequence
    requires WritableAlphabet<typename t::seq_alphabet>;
    requires WritableAlphabet<typename t::seq_legal_alphabet>;
    requires ExplicitlyConvertibleTo<typename t::seq_legal_alphabet, typename t::seq_alphabet>;
    requires SequenceContainer<typename t::template seq_container<typename t::seq_alphabet>>;
//    requires SequenceContainer
//        <typename t::template seq_container_container
//            <typename t::template seq_container
//                <typename t::seq_alphabet>>>;

    // id
    requires WritableAlphabet<typename t::id_alphabet>;
    requires SequenceContainer<typename t::template id_container<typename t::id_alphabet>>;
//    requires SequenceContainer
//        <typename t::template id_container_container
//            <typename t::template id_container
//                <typename t::id_alphabet>>>;

    // bpp
    requires std::is_floating_point_v<typename t::bpp_prob>;
    requires std::numeric_limits<typename t::bpp_partner>::is_integer;

//    requires Container // TODO check Associative Container Concept when implemented
//        <typename t::template bpp_queue
//            <typename t::template bpp_item
//                <typename t::bpp_prob, typename t::bpp_partner>>>
//        && requires(typename t::template bpp_queue // TODO maybe implement also a version that allows emplace_back
//            <typename t::template bpp_item
//                <typename t::bpp_prob, typename t::bpp_partner>> value) { value.emplace(1.0, 1); };
//    requires SequenceContainer
//        <typename t::template bpp_container
//            <typename t::template bpp_queue
//                 <typename t::template bpp_item
//                      <typename t::bpp_prob, typename t::bpp_partner>>>>;
//    requires SequenceContainer
//        <typename t::template bpp_container_container
//            <typename t::template bpp_container
//                <typename t::template bpp_queue
//                    <typename t::template bpp_item
//                        <typename t::bpp_prob, typename t::bpp_partner>>>>>;

    // structure
    requires std::is_same_v<typename t::structure_alphabet, dssp9> // TODO(joergi-w) add aa_structure_concept
          || RnaStructureAlphabet<typename t::structure_alphabet>;
    requires SequenceContainer<typename t::template structure_container<typename t::structure_alphabet>>;
//    requires SequenceContainer
//        <typename t::template structure_container_container
//            <typename t::template structure_container
//                <typename t::structure_alphabet>>>;

    // structured sequence: tuple composites of seq and structure
    requires std::is_base_of_v<alphabet_tuple_base
        <typename t::template structured_seq_alphabet
            <typename t::seq_alphabet, typename t::structure_alphabet>,
             typename t::seq_alphabet, typename t::structure_alphabet>,
        typename t::template structured_seq_alphabet<typename t::seq_alphabet, typename t::structure_alphabet>>;
//    requires SequenceContainer
//        <typename t::template structured_seq_container
//            <typename t::template structured_seq_alphabet
//                <typename t::seq_alphabet, typename t::structure_alphabet>>>;
//    requires SequenceContainer
//        <typename t::template structured_seq_container_container
//            <typename t::template structured_seq_container
//                <typename t::template structured_seq_alphabet
//                    <typename t::seq_alphabet, typename t::structure_alphabet>>>>;

    // energy: std::optional of floating point number
    requires std::is_floating_point_v<typename t::energy_type::value_type>;
    requires SequenceContainer<typename t::template energy_container<typename t::energy_type>>;

    // reactivity [error]
    requires std::is_floating_point_v<typename t::react_type>;
    requires SequenceContainer<typename t::template react_container<typename t::react_type>>;
//    requires SequenceContainer
//        <typename t::template react_container_container
//            <typename t::template react_container
//                <typename t::react_type>>>;

    // comment
    requires WritableAlphabet<typename t::comment_alphabet>;
    requires SequenceContainer<typename t::template comment_container<typename t::comment_alphabet>>;
//    requires SequenceContainer
//        <typename t::template comment_container_container
//            <typename t::template comment_container
//                <typename t::comment_alphabet>>>;

    // offset
    requires std::numeric_limits<typename t::offset_type>::is_integer;
    requires SequenceContainer<typename t::template offset_container<typename t::offset_type>>;
};
//!\endcond

// ----------------------------------------------------------------------------
// structure_file_input_default_traits
// ----------------------------------------------------------------------------

/*!\brief The default traits for seqan3::structure_file_input
 * \implements StructureFileInputTraits
 * \ingroup structure_file
 *
 * \details
 *
 * If you wish to change a single or a few types from the default, just inherit from this class and
 * "overwrite" the respective type definitions.
 *
 * This example will make the file read into a smaller alphabet:
 *
 * \snippet test/unit/io/structure_file/structure_file_input_test.cpp structure_file_input_class mod_traits
 */
struct structure_file_input_default_traits_rna
{
    /*!\name Member types
     * \brief Definitions to satisfy seqan3::StructureFileInputTraits.
     * \{
     */

    // sequence

    //!\brief The sequence alphabet is seqan3::rna5.
    using seq_alphabet                       = rna5;

    //!\brief The legal sequence alphabet for parsing is seqan3::rna15.
    using seq_legal_alphabet                 = rna15;

    //!\brief The type of an RNA sequence is std::vector.
    template<typename _seq_alphabet>
    using seq_container                      = std::vector<_seq_alphabet>;

    //!\brief The container for sequences is seqan3::concatenated_sequences.
    template<typename _seq_container>
    using seq_container_container            = concatenated_sequences<_seq_container>;

    // id

    //!\brief The alphabet for an identifier string is char.
    using id_alphabet                        = char;

    //!\brief The string type for an identifier is std::basic_string.
    template<typename _id_alphabet>
    using id_container                       = std::basic_string<_id_alphabet>;

    //!\brief The container for identifier strings is seqan3::concatenated_sequences.
    template<typename _id_container>
    using id_container_container             = concatenated_sequences<_id_container>;

    // base pair probability structure

    //!\brief The type for a base pair probability is double.
    using bpp_prob                           = double;

    //!\brief The type for the partner position of a base pair probability is size_t.
    using bpp_partner                        = size_t;

    //!\brief The type of a base pair item is std::pair<double, size_t>.
    template<typename _bpp_prob, typename _bpp_partner>
    using bpp_item                           = std::pair<_bpp_prob, _bpp_partner>;

    //!\brief A queue of base pair items sorted by probability is realised with std::set.
    template<typename _bpp_item>
    using bpp_queue                          = std::set<_bpp_item>;

    //!\brief A string over all bases containing the respective interaction queues is represented as std::vector.
    template<typename _bpp_queue>
    using bpp_container                      = std::vector<_bpp_queue>;

    //!\brief The container for interaction strings is std::vector.
    template<typename _bpp_container>
    using bpp_container_container            = std::vector<_bpp_container>;

    // fixed structure

    //!\brief The alphabet for a structure annotation is seqan3::phred42.
    using structure_alphabet                 = wuss51;

    //!\brief The string type for a structure annotation is std::vector.
    template<typename _structure_alphabet>
    using structure_container                = std::vector<_structure_alphabet>;

    //!\brief The container for structure annotation strings is seqan3::concatenated_sequences.
    template<typename _structure_container>
    using structure_container_container      = concatenated_sequences<_structure_container>;

    // combined sequence and structure

    //!\brief The combined structured sequence alphabet is seqan3::structured_rna<seqan3::rna5, seqan3::wuss51>.
    template<typename _seq_alphabet, typename _structure_alphabet>
    using structured_seq_alphabet            = structured_rna<_seq_alphabet, _structure_alphabet>;

    //!\brief The type of a structured RNA sequence is std::vector.
    template<typename _structured_seq_alphabet>
    using structured_seq_container           = std::vector<_structured_seq_alphabet>;

    //!\brief The container for sequences is seqan3::concatenated_sequences.
    template<typename _structured_seq_container>
    using structured_seq_container_container = concatenated_sequences<_structured_seq_container>;

    // energy

    //!\brief The type of the energy is std::optional<double>.
    using energy_type                        = std::optional<double>;

    //!\brief The type of a container of energy values is std::vector.
    template<typename _energy_type>
    using energy_container                   = std::vector<_energy_type>;

    // reactivity [error]

    //!\brief The type of the reactivity and reactivity error is double.
    using react_type                         = double;

    //!\brief The type of a string of reactivity values is std::vector.
    template<typename _react_type>
    using react_container                    = std::vector<_react_type>;

    //!\brief The type of a container of reactivity strings is std::vector.
    template<typename _react_container>
    using react_container_container          = std::vector<_react_container>;

    // comment

    //!\brief The alphabet for a comment string is char.
    using comment_alphabet                   = char;

    //!\brief The string type for a comment is std::basic_string.
    template<typename _comment_alphabet>
    using comment_container                  = std::basic_string<_comment_alphabet>;

    //!\brief The container for comments is seqan3::concatenated_sequences.
    template<typename _comment_container>
    using comment_container_container        = concatenated_sequences<_comment_container>;

    // offset

    //!\brief The type of the offset is size_t.
    using offset_type                        = size_t;

    //!\brief The type of a container of offset values is std::vector.
    template<typename _offset_type>
    using offset_container                   = std::vector<_offset_type>;
    //!\}
};

//!\brief A traits type that specifies input as amino acids.
//!\ingroup structure_file
struct structure_file_input_default_traits_aa : structure_file_input_default_traits_rna
{
    /*!\name Member types
     * \brief Definitions to satisfy seqan3::StructureFileInputTraits.
     * \{
     */

    //!\brief The sequence alphabet is seqan3::aa27.
    using seq_alphabet             = aa27;
    //!\brief The legal sequence alphabet for parsing is seqan3::aa27.
    using seq_legal_alphabet       = aa27;
    //!\brief The structure annotation alphabet is seqan3::dssp9.
    using structure_alphabet       = dssp9;
    //!\brief The combined structured sequence alphabet is seqan3::structured_aa<seqan3::aa27, seqan3::dssp9>.
    template<typename _seq_alphabet, typename _structure_alphabet>
    using structured_seq_alphabet  = structured_aa<_seq_alphabet, _structure_alphabet>;
    //!\}
};

// ----------------------------------------------------------------------------
// structure_file_input
// ----------------------------------------------------------------------------

/*!\brief A class for reading structured sequence files, e.g. Stockholm, Connect, Vienna, ViennaRNA bpp matrix ...
 * \ingroup structure_file
 * \tparam traits_type        An auxiliary type that defines certain member types and constants, must satisfy
 *                            seqan3::StructureFileInputTraits.
 * \tparam selected_field_ids A seqan3::fields type with the list and order of desired record entries; all fields
 *                            must be in seqan3::structure_file_input::field_ids.
 * \tparam valid_formats      A seqan3::type_list of the selectable formats (each must meet
 *                            seqan3::StructureFileInputFormat).
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
 * The structured sequence file abstraction supports reading ten different fields:
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
 * The first three fields are retrieved by default (and in that order). The seqan3::field::STRUCTURED_SEQ may be
 * selected to have sequence and structure directly stored in a more memory-efficient combined container.
 * If you select this field you must not select seqan3::field::SEQ or seqan3::field::STRUCTURE.
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
 * \snippet test/snippet/io/structure_file/structure_file_input.cpp auto_temp_deduc
 *
 * Reading from an std::istringstream:
 *
 * \snippet test/snippet/io/structure_file/structure_file_input.cpp stringstream_read
 *
 * Note that this is not the same as writing `structure_file_input<>` (with angle brackets). In the latter case they are
 * explicitly set to their default values, in the former case
 * [automatic deduction](http://en.cppreference.com/w/cpp/language/class_template_argument_deduction) happens which
 * chooses different parameters depending on the constructor arguments. For opening from file, `structure_file_input<>`
 * would have also worked, but for opening from stream it would not have.
 *
 * In some cases, you do need to specify the arguments, e.g. if you want to read amino acids:
 *
 * \snippet test/snippet/io/structure_file/structure_file_input.cpp arg_spec
 *
 * You can define your own traits type to further customise the types used by and returned by this class, see
 * seqan3::structure_file_default_traits_rna for more details. As mentioned above, specifying at least one
 * template parameter yourself means that you loose automatic deduction so if you want to read amino acids **and**
 * want to read from a string stream you need to give all types yourself:
 *
 * \snippet test/snippet/io/structure_file/structure_file_input.cpp trait_def
 *
 * ### Reading record-wise
 *
 * You can iterate over this file record-wise:
 *
 * \snippet test/snippet/io/structure_file/structure_file_input.cpp record_iter
 *
 * In the above example, rec has the type \ref record_type which is a specialisation of seqan3::record and behaves
 * like an std::tuple (that's why we can access it via get). Instead of using the seqan3::field based interface on
 * the record, you could also use `std::get<0>` or even `std::get<rna4_vector>` to retrieve the sequence, but it is
 * not recommended, because it is more error-prone.
 *
 * *Note:* It is important to write `auto &` and not just `auto`, otherwise you will copy the record on every iteration.
 * Since the buffer gets "refilled" on every iteration, you can also move the data out of the record if you want
 * to store it somewhere without copying:
 *
 * \snippet test/snippet/io/structure_file/structure_file_input.cpp data_out
 *
 * ### Reading record-wise (decomposed records)
 *
 * Instead of using `get` on the record, you can also use
 * [structured bindings](http://en.cppreference.com/w/cpp/language/structured_binding)
 * to decompose the record into its elements:
 *
 * \snippet test/snippet/io/structure_file/structure_file_input.cpp structured_bindings
 *
 * In this case you immediately get the three elements of the tuple: `seq` of \ref seq_type, `id` of
 * \ref id_type and `structure` of \ref structure_type. **But beware: with structured bindings you do need
 * to get the order of elements correctly!**
 *
 * ### Reading record-wise (custom fields)
 *
 * If you want to skip specific fields from the record you can pass a non-empty fields trait object to the
 * structure_file_input constructor to select the fields that should be read from the input. For example to choose a
 * combined field for SEQ and STRUCTURE (see above). Or to never actually read the STRUCTURE, if you don't need it.
 * The following snippets demonstrate the usage of such a fields trait object.
 *
 * \snippet test/snippet/io/structure_file/structure_file_input.cpp skip_fields
 *
 * When reading a file, all fields not present in the file (but requested implicitly or via the `selected_field_ids`
 * parameter) are ignored.
 *
 * ### Views on files
 *
 * Since SeqAn files are ranges, you can also create views over files. A useful example is to filter the records
 * based on certain criteria, e.g. minimum length of the sequence field:
 *
 * \snippet test/snippet/io/structure_file/structure_file_input.cpp filter_criteria
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
 * \snippet test/snippet/io/structure_file/structure_file_input.cpp data_storage
 * \snippet test/snippet/io/structure_file/structure_file_input.cpp col_read
 *
 * Note that for this to make sense, your storage data types need to be identical to the corresponding column types
 * of the file. If you require different column types you can specify you own traits, see
 * seqan3::StructureFileInputTraits.
 *
 * ### Formats
 *
 * Currently, the only implemented format is seqan3::format_vienna. More formats will follow soon.
 */
template<StructureFileInputTraits traits_type_ = structure_file_input_default_traits_rna,
         detail::Fields selected_field_ids_ = fields<field::SEQ, field::ID, field::STRUCTURE>,
         detail::TypeListOfStructureFileInputFormats valid_formats_
             = type_list<format_vienna>,
         Char stream_char_type_ = char>
class structure_file_input
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
    //!\brief Character type of the stream(s), usually `char`.
    using stream_char_type   = stream_char_type_;
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
                  "of structure_file_input::field_ids for the accepted values.");

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
                <typename traits_type::bpp_prob, typename traits_type::bpp_partner>>>;
    //!\brief The type of the structure field (default std::vector of seqan3::wuss51).
    using structure_type      = typename traits_type::template structure_container
        <typename traits_type::structure_alphabet>;
    //!\brief The type of the sequence-structure field (default std::vector of structured_rna<rna5, wuss51>).
    using structured_seq_type = typename traits_type::template structured_seq_container
        <typename traits_type::template structured_seq_alphabet
            <typename traits_type::seq_alphabet, typename traits_type::structure_alphabet>>;
    //!\brief The type of the energy field (default double).
    using energy_type         = typename traits_type::energy_type;
    //!\brief The type of the reactivity and reactivity error fields (default double).
    using react_type          = typename traits_type::template react_container<typename traits_type::react_type>;
    //!\brief The type of the comment field (default double).
    using comment_type        = typename traits_type::template comment_container
        <typename traits_type::comment_alphabet>;
    //!\brief The type of the offset field (default size_t).
    using offset_type         = typename traits_type::offset_type;

    //!\brief The previously defined types aggregated in a seqan3::type_list.
    using field_types    = type_list<seq_type, id_type, bpp_type, structure_type, structured_seq_type, energy_type,
                                     react_type, react_type, comment_type, offset_type>;

    //!\brief The type of the record, a specialisation of seqan3::record; acts as a tuple of the selected field types.
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
    using iterator        = detail::in_file_iterator<structure_file_input>;
    //!\brief The const iterator type is void, because files are not const-iterable.
    using const_iterator  = void;
    //!\brief The type returned by end().
    using sentinel        = std::ranges::default_sentinel_t;
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    //!\brief Default constructor is explicitly deleted, you need to give a stream or file name.
    structure_file_input()                                         = delete;
    //!\brief Copy construction is explicitly deleted, because you cannot have multiple access to the same file.
    structure_file_input(structure_file_input const &)             = delete;
    //!\brief Copy assignment is explicitly deleted, because you cannot have multiple access to the same file.
    structure_file_input & operator=(structure_file_input const &) = delete;
    //!\brief Move construction is defaulted.
    structure_file_input(structure_file_input &&)                  = default;
    //!\brief Move assignment is defaulted.
    structure_file_input & operator=(structure_file_input &&)      = default;
    //!\brief Destructor is defaulted.
    ~structure_file_input()                                        = default;

    /*!\brief Construct from filename.
     * \param[in] filename Path to the file you wish to open.
     * \param[in] fields_tag A seqan3::fields tag. [optional]
     * \throws seqan3::file_open_error if the file could not be opened, e.g. non-existent, non-readable, unknown format.
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
    structure_file_input(std::filesystem::path filename,
                         selected_field_ids const & SEQAN3_DOXYGEN_ONLY(fields_tag) = selected_field_ids{}) :
        primary_stream{new std::ifstream{filename, std::ios_base::in | std::ios::binary}, stream_deleter_default}
    {
        if (!primary_stream->good())
            throw file_open_error{"Could not open file " + filename.string() + " for reading."};

        // possibly add intermediate decompression stream
        secondary_stream = detail::make_secondary_istream(*primary_stream, filename);

        // initialise format handler
        detail::set_format(format, filename);

        // buffer first record
        read_next_record();
    }

    /*!\brief Construct from an existing stream and with specified format.
     * \tparam file_format The format of the file in the stream, must satisfy
     * seqan3::StructureFileInputFormat.
     * \param[in] stream The stream to operate on; must be derived of std::basic_istream.
     * \param[in] format_tag The file format tag.
     * \param[in] fields_tag A seqan3::fields tag. [optional]
     *
     * \details
     *
     * ### Decompression
     *
     * This constructor transparently applies a decompression stream on top of the stream in case
     * the file is detected as being compressed.
     * See the section on \link io_compression compression and decompression \endlink for more information.
     */
    template<IStream2 stream_t, StructureFileInputFormat file_format>
    structure_file_input(stream_t & stream,
                         file_format const & SEQAN3_DOXYGEN_ONLY(format_tag),
                         selected_field_ids const & SEQAN3_DOXYGEN_ONLY(fields_tag) = selected_field_ids{}) :
        primary_stream{&stream, stream_deleter_noop},
        format{detail::structure_file_input_format<file_format>{}}
    {
        static_assert(meta::in<valid_formats, file_format>::value,
                      "You selected a format that is not in the valid_formats of this file.");

        // possibly add intermediate decompression stream
        secondary_stream = detail::make_secondary_istream(*primary_stream);

        // buffer first record
        read_next_record();
    }

    //!\overload
    template<IStream2 stream_t, StructureFileInputFormat file_format>
    structure_file_input(stream_t && stream,
                         file_format const & SEQAN3_DOXYGEN_ONLY(format_tag),
                         selected_field_ids const & SEQAN3_DOXYGEN_ONLY(fields_tag) = selected_field_ids{}) :
        primary_stream{new stream_t{std::move(stream)}, stream_deleter_default},
        format{detail::structure_file_input_format<file_format>{}}
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
     * \snippet test/snippet/io/structure_file/structure_file_input.cpp ref_return
     *
     * It most situations using the iterator interface or a range-based for-loop are preferable to using front(),
     * because you can only move to the next record via the iterator.
     *
     * In any case, don't forget the reference! If you want to save the data from the record elsewhere, use move:
     *
     * \snippet test/snippet/io/structure_file/structure_file_input.cpp move
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
    friend auto & get(structure_file_input & file)
    {
        static_assert(structure_file_input::selected_field_ids::contains(f),
                      "You requested a field via get that was not selected for the file.");

        file.read_columns();

        return seqan3::get<f>(file.columns_buffer);
    }

    //!\copydoc get
    template<field f>
    friend auto && get(structure_file_input && file)
    {
        return std::move(get<f>(file));
    }

    //!\copydoc get
    template<size_t i>
    friend auto & get(structure_file_input & file)
    {
        static_assert(i < structure_file_input::selected_field_ids::as_array.size(),
                      "You requested a field number larger than the number of selected fields for the file.");
        file.read_columns();

        return std::get<i>(file.columns_buffer);
    }

    //!\copydoc get
    template<size_t i>
    friend auto && get(structure_file_input && file)
    {
        return std::move(get<i>(file));
    }

    //!\copydoc get
    template<typename t>
    friend auto & get(structure_file_input & file)
    {
        file.read_columns();

        return std::get<t>(file.columns_buffer);
    }

    //!\copydoc get
    template<typename t>
    friend auto && get(structure_file_input && file)
    {
        return std::move(get<t>(file));
    }
    //!\}

    //!\brief The options are public and its members can be set directly.
    structure_file_input_options<typename traits_type::seq_legal_alphabet,
                                 selected_field_ids::contains(field::STRUCTURED_SEQ)> options;

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
    using format_type = typename detail::variant_from_tags<valid_formats, detail::structure_file_input_format>::type;
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
            if constexpr (selected_field_ids::contains(field::STRUCTURED_SEQ))
            {
                static_assert(!selected_field_ids::contains(field::STRUCTURE),
                              "You may not select field::STRUCTURED_SEQ and field::STRUCTURE at the same time.");
                static_assert(!selected_field_ids::contains(field::SEQ),
                              "You may not select field::STRUCTURED_SEQ and field::SEQ at the same time.");
                f.read(*secondary_stream,
                       options,
                       detail::get_or_ignore<field::STRUCTURED_SEQ>(record_buffer), // seq
                       detail::get_or_ignore<field::ID>(record_buffer),
                       detail::get_or_ignore<field::BPP>(record_buffer),
                       detail::get_or_ignore<field::STRUCTURED_SEQ>(record_buffer), // structure
                       detail::get_or_ignore<field::ENERGY>(record_buffer),
                       detail::get_or_ignore<field::REACT>(record_buffer),
                       detail::get_or_ignore<field::REACT_ERR>(record_buffer),
                       detail::get_or_ignore<field::COMMENT>(record_buffer),
                       detail::get_or_ignore<field::OFFSET>(record_buffer));
            }
            else
            {
                f.read(*secondary_stream,
                       options,
                       detail::get_or_ignore<field::SEQ>(record_buffer),
                       detail::get_or_ignore<field::ID>(record_buffer),
                       detail::get_or_ignore<field::BPP>(record_buffer),
                       detail::get_or_ignore<field::STRUCTURE>(record_buffer),
                       detail::get_or_ignore<field::ENERGY>(record_buffer),
                       detail::get_or_ignore<field::REACT>(record_buffer),
                       detail::get_or_ignore<field::REACT_ERR>(record_buffer),
                       detail::get_or_ignore<field::COMMENT>(record_buffer),
                       detail::get_or_ignore<field::OFFSET>(record_buffer));
            }
        }, format);
    }

    //!\brief Read the entire file into the internal column buffers.
    void read_columns()
    {
        //TODO don't do multiple visits
        //TODO create specialised version for concatenated_sequences where we append on the concat
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

    //!\brief Befriend iterator so it can access the buffers.
    friend iterator;
};

/*!\name Type deduction guides
 * \relates seqan3::structure_file_input
 * \{
 */

//!\brief Deduction of the selected fields, the file format and the stream type.
template <IStream2                 stream_type,
          StructureFileInputFormat file_format,
          detail::Fields           selected_field_ids>
structure_file_input(stream_type && stream, file_format const &, selected_field_ids const &)
    -> structure_file_input<typename structure_file_input<>::traits_type,       // actually use the default
                            selected_field_ids,
                            type_list<file_format>,
                            typename std::remove_reference_t<stream_type>::char_type>;

//!\overload
template <IStream2                 stream_type,
          StructureFileInputFormat file_format,
          detail::Fields           selected_field_ids>
structure_file_input(stream_type & stream, file_format const &, selected_field_ids const &)
    -> structure_file_input<typename structure_file_input<>::traits_type,       // actually use the default
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
 * \ingroup structure_file
 * \see std::tuple_size_v
 */
template<seqan3::StructureFileInputTraits                    traits_type,
         seqan3::detail::Fields                              selected_field_ids,
         seqan3::detail::TypeListOfStructureFileInputFormats valid_formats,
         seqan3::Char                                        stream_char_t>
struct tuple_size<seqan3::structure_file_input<traits_type, selected_field_ids, valid_formats, stream_char_t>>
{
    //!\brief The value equals the number of selected fields in the file.
    static constexpr size_t value = selected_field_ids::as_array.size();
};

/*!\brief Obtains the type of the specified element.
 * \implements seqan3::TransformationTrait
 * \ingroup structure_file
 * \see [std::tuple_element](https://en.cppreference.com/w/cpp/utility/tuple/tuple_element)
 */
template<size_t                                              elem_no,
         seqan3::StructureFileInputTraits                    traits_type,
         seqan3::detail::Fields                              selected_field_ids,
         seqan3::detail::TypeListOfStructureFileInputFormats valid_formats,
         seqan3::Char                                        stream_char_t>
struct tuple_element<elem_no,
                     seqan3::structure_file_input<traits_type, selected_field_ids, valid_formats, stream_char_t>>
    : tuple_element<elem_no, typename seqan3::structure_file_input<traits_type,
                                                                   selected_field_ids,
                                                                   valid_formats,
                                                                   stream_char_t>::file_as_tuple_type>
{};

} // namespace std

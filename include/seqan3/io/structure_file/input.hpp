// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

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
#include <seqan3/alphabet/aminoacid/all.hpp>
#include <seqan3/alphabet/nucleotide/rna15.hpp>
#include <seqan3/alphabet/nucleotide/rna5.hpp>
#include <seqan3/alphabet/structure/all.hpp>
#include <seqan3/core/metafunction/basic.hpp>
#include <seqan3/io/stream/concept.hpp>
#include <seqan3/io/exception.hpp>
#include <seqan3/io/filesystem.hpp>
#include <seqan3/io/record.hpp>
#include <seqan3/io/detail/in_file_iterator.hpp>
#include <seqan3/io/detail/record.hpp>
#include <seqan3/io/structure_file/input_format_concept.hpp>
#include <seqan3/io/structure_file/input_options.hpp>
#include <seqan3/io/structure_file/format_vienna.hpp>
#include <seqan3/range/container/concatenated_sequences.hpp>

namespace seqan3
{
// ----------------------------------------------------------------------------
// structure_file_input_traits_concept
// ----------------------------------------------------------------------------

/*!\interface seqan3::structure_file_input_traits_concept <>
 * \brief The requirements a traits_type for seqan3::structure_file_in must meet.
 * \ingroup structure_file
 */
/*!\name Requirements for seqan3::structure_file_input_traits_concept
 * \brief You can expect these **member types** of all types that satisfy seqan3::structure_file_input_traits_concept.
 * \memberof seqan3::structure_file_input_traits_concept
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
 * \memberof seqan3::structure_file_input_traits_concept
 * \brief Alphabet of the characters for the seqan3::field::SEQ; must satisfy seqan3::alphabet_concept.
 */
/*!\typedef using seq_legal_alphabet
 * \memberof seqan3::structure_file_input_traits_concept
 * \brief Intermediate alphabet for seqan3::field::SEQ; must satisfy seqan3::alphabet_concept and be convertible to
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
 * \memberof seqan3::structure_file_input_traits_concept
 * \brief Type template of the seqan3::field::SEQ, a container template over `seq_alphabet`;
 * must satisfy seqan3::sequence_container_concept.
 */
/*!\typedef using seq_container_container
 * \memberof seqan3::structure_file_input_traits_concept
 * \brief Type template of a column of seqan3::field::SEQ, a container template that can hold multiple
 * `seq_container`; must satisfy seqan3::sequence_container_concept.
 */
/*!\typedef using id_alphabet
 * \memberof seqan3::structure_file_input_traits_concept
 * \brief Alphabet of the characters for the seqan3::field::ID; must satisfy seqan3::alphabet_concept.
 */
/*!\typedef using id_container
 * \memberof seqan3::structure_file_input_traits_concept
 * \brief Type template of the seqan3::field::ID, a container template over `id_alphabet`;
 * must satisfy seqan3::sequence_container_concept.
 */
/*!\typedef using id_container_container
 * \memberof seqan3::structure_file_input_traits_concept
 * \brief Type template of a column of seqan3::field::ID, a container template that can hold multiple
 * `id_container`; must satisfy seqan3::sequence_container_concept.
 */
/*!\typedef using bpp_prob
 * \memberof seqan3::structure_file_input_traits_concept
 * \brief Data type for the base pair probabilities in seqan3::field::BPP; must satisfy std::is_floating_point.
 */
/*!\typedef using bpp_partner
 * \memberof seqan3::structure_file_input_traits_concept
 * \brief Data type for the partner index of an interaction in seqan3::field::BPP; must satisfy
 * std::numeric_limits::is_integer.
 */
/*!\typedef using bpp_queue
 * \memberof seqan3::structure_file_input_traits_concept
 * \brief A container template representing a set of interactions of type bpp_item,
 * which are (comparable) tuples of `bpp_prob` and `bpp_partner`;
 * must satisfy seqan3::container_concept and must provide an std::emplace(bpp_prob, bpp_partner) function.
 */
/*!\typedef using bpp_container
 * \memberof seqan3::structure_file_input_traits_concept
 * \brief Type template of the seqan3::field::BPP, a container template over a set (bpp_queue) of interactions;
 * must satisfy seqan3::sequence_container_concept.
 */
/*!\typedef using bpp_container_container
 * \memberof seqan3::structure_file_input_traits_concept
 * \brief Type template of a column of seqan3::field::BPP, a container template that can hold multiple
 * `bpp_container`; must satisfy seqan3::sequence_container_concept.
 */
/*!\typedef using structure_alphabet
 * \memberof seqan3::structure_file_input_traits_concept
 * \brief Alphabet of the characters for the seqan3::field::STRUCTURE; must satisfy seqan3::rna_structure_concept.
 */
/*!\typedef using structure_container
 * \memberof seqan3::structure_file_input_traits_concept
 * \brief Type template of the seqan3::field::STRUCTURE, a container template over `structure_alphabet`;
 * must satisfy seqan3::sequence_container_concept.
 */
/*!\typedef using structure_container_container
 * \memberof seqan3::structure_file_input_traits_concept
 * \brief Type template of a column of seqan3::field::STRUCTURE, a container template that can hold multiple
 * `structure_container`; must satisfy seqan3::sequence_container_concept.
 */
/*!\typedef using energy_type
 * \memberof seqan3::structure_file_input_traits_concept
 * \brief Type template of the seqan3::field::ENERGY; must be std::optional of a type satisfying std::is_floating_point.
 * \details
 * If the file record contains an energy, the value can be obtained through dereference or std::optional::value.
 * Otherwise, operator bool or std::optional::has_value return false.
 */
/*!\typedef using energy_container
 * \memberof seqan3::structure_file_input_traits_concept
 * \brief Type template of a column of seqan3::field::ENERGY, a container template that can hold multiple `energy_type`;
 * must satisfy seqan3::sequence_container_concept.
 */
/*!\typedef using react_type
 * \memberof seqan3::structure_file_input_traits_concept
 * \brief Data type for the reactivity and reactivity error in seqan3::field::REACT and seqan3::field::REACT_ERR,
 * respectively; must satisfy std::is_floating_point.
 */
/*!\typedef using react_container
 * \memberof seqan3::structure_file_input_traits_concept
 * \brief Type template of the seqan3::field::REACT and seqan3::field::REACT_ERR, a container template over
 * `react_type`; must satisfy seqan3::sequence_container_concept.
 */
/*!\typedef using react_container_container
 * \memberof seqan3::structure_file_input_traits_concept
 * \brief Type template of a column of seqan3::field::REACT and seqan3::field::REACT_ERR, a container template that
 * can hold multiple `react_container`; must satisfy seqan3::sequence_container_concept.
 */
/*!\typedef using comment_alphabet
 * \memberof seqan3::structure_file_input_traits_concept
 * \brief Alphabet of the characters for the seqan3::field::COMMENT; must satisfy seqan3::alphabet_concept.
 */
/*!\typedef using comment_container
 * \memberof seqan3::structure_file_input_traits_concept
 * \brief Type template of the seqan3::field::COMMENT, a container template over `comment_alphabet`;
 * must satisfy seqan3::sequence_container_concept.
 */
/*!\typedef using comment_container_container
 * \memberof seqan3::structure_file_input_traits_concept
 * \brief Type template of a column of seqan3::field::COMMENT, a container template that can hold multiple
 * `comment_container`; must satisfy seqan3::sequence_container_concept.
 */
/*!\typedef using offset_type
 * \memberof seqan3::structure_file_input_traits_concept
 * \brief Type template of the seqan3::field::OFFSET; must statisfy std::numeric_limits::is_integer.
 */
/*!\typedef using offset_container
 * \memberof seqan3::structure_file_input_traits_concept
 * \brief Type template of a column of seqan3::field::OFFSET, a container template that can hold multiple `offset_type`;
 * must satisfy seqan3::sequence_container_concept.
 */
//!\}
//!\cond
template<typename t>
SEQAN3_CONCEPT structure_file_input_traits_concept = requires(t v)
{
    // TODO(joergi-w) The expensive concept checks are currently omitted. Check again when compiler has improved.
    // sequence
    requires alphabet_concept<typename t::seq_alphabet>;
    requires alphabet_concept<typename t::seq_legal_alphabet>;
    requires explicitly_convertible_to_concept<typename t::seq_legal_alphabet, typename t::seq_alphabet>;
    requires sequence_container_concept<typename t::template seq_container<typename t::seq_alphabet>>;
//    requires sequence_container_concept
//        <typename t::template seq_container_container
//            <typename t::template seq_container
//                <typename t::seq_alphabet>>>;

    // id
    requires alphabet_concept<typename t::id_alphabet>;
    requires sequence_container_concept<typename t::template id_container<typename t::id_alphabet>>;
//    requires sequence_container_concept
//        <typename t::template id_container_container
//            <typename t::template id_container
//                <typename t::id_alphabet>>>;

    // bpp
    requires std::is_floating_point_v<typename t::bpp_prob>;
    requires std::numeric_limits<typename t::bpp_partner>::is_integer;

//    requires container_concept // TODO check Associative Container Concept when implemented
//        <typename t::template bpp_queue
//            <typename t::template bpp_item
//                <typename t::bpp_prob, typename t::bpp_partner>>>
//        && requires(typename t::template bpp_queue // TODO maybe implement also a version that allows emplace_back
//            <typename t::template bpp_item
//                <typename t::bpp_prob, typename t::bpp_partner>> value) { value.emplace(1.0, 1); };
//    requires sequence_container_concept
//        <typename t::template bpp_container
//            <typename t::template bpp_queue
//                 <typename t::template bpp_item
//                      <typename t::bpp_prob, typename t::bpp_partner>>>>;
//    requires sequence_container_concept
//        <typename t::template bpp_container_container
//            <typename t::template bpp_container
//                <typename t::template bpp_queue
//                    <typename t::template bpp_item
//                        <typename t::bpp_prob, typename t::bpp_partner>>>>>;

    // structure
    requires std::is_same_v<typename t::structure_alphabet, dssp9> // TODO(joergi-w) add aa_structure_concept
          || rna_structure_concept<typename t::structure_alphabet>;
    requires sequence_container_concept<typename t::template structure_container<typename t::structure_alphabet>>;
//    requires sequence_container_concept
//        <typename t::template structure_container_container
//            <typename t::template structure_container
//                <typename t::structure_alphabet>>>;

    // structured sequence: cartesian compositions of seq and structure
    requires std::is_base_of_v<cartesian_composition
        <typename t::template structured_seq_alphabet
            <typename t::seq_alphabet, typename t::structure_alphabet>,
             typename t::seq_alphabet, typename t::structure_alphabet>,
        typename t::template structured_seq_alphabet<typename t::seq_alphabet, typename t::structure_alphabet>>;
//    requires sequence_container_concept
//        <typename t::template structured_seq_container
//            <typename t::template structured_seq_alphabet
//                <typename t::seq_alphabet, typename t::structure_alphabet>>>;
//    requires sequence_container_concept
//        <typename t::template structured_seq_container_container
//            <typename t::template structured_seq_container
//                <typename t::template structured_seq_alphabet
//                    <typename t::seq_alphabet, typename t::structure_alphabet>>>>;

    // energy: std::optional of floating point number
    requires std::is_floating_point_v<typename t::energy_type::value_type>;
    requires sequence_container_concept<typename t::template energy_container<typename t::energy_type>>;

    // reactivity [error]
    requires std::is_floating_point_v<typename t::react_type>;
    requires sequence_container_concept<typename t::template react_container<typename t::react_type>>;
//    requires sequence_container_concept
//        <typename t::template react_container_container
//            <typename t::template react_container
//                <typename t::react_type>>>;

    // comment
    requires alphabet_concept<typename t::comment_alphabet>;
    requires sequence_container_concept<typename t::template comment_container<typename t::comment_alphabet>>;
//    requires sequence_container_concept
//        <typename t::template comment_container_container
//            <typename t::template comment_container
//                <typename t::comment_alphabet>>>;

    // offset
    requires std::numeric_limits<typename t::offset_type>::is_integer;
    requires sequence_container_concept<typename t::template offset_container<typename t::offset_type>>;
};
//!\endcond

// ----------------------------------------------------------------------------
// structure_file_input_default_traits
// ----------------------------------------------------------------------------

/*!\brief The default traits for seqan3::structure_file_in
 * \implements structure_file_input_traits_concept
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
     * \brief Definitions to satisfy seqan3::structure_file_input_traits_concept.
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
    using bpp_prob                           = double;
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
    template<typename _structure_alphabet>
    using structure_container                = std::vector<_structure_alphabet>;
    template<typename _structure_container>
    using structure_container_container      = concatenated_sequences<_structure_container>;

    // combined sequence and structure
    template<typename _seq_alphabet, typename _structure_alphabet>
    using structured_seq_alphabet            = structured_rna<_seq_alphabet, _structure_alphabet>;
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

//!\brief A traits type that specifies input as amino acids.
//!\ingroup structure_file
struct structure_file_input_default_traits_aa : structure_file_input_default_traits_rna
{
    /*!\name Member types
     * \brief Definitions to satisfy seqan3::structure_file_input_traits_concept.
     * \{
     */
    using seq_alphabet             = aa27;
    using seq_legal_alphabet       = aa27;
    using structure_alphabet       = dssp9;
    template<typename _seq_alphabet, typename _structure_alphabet>
    using structured_seq_alphabet  = structured_aa<_seq_alphabet, _structure_alphabet>;
    //!\}
};

// ----------------------------------------------------------------------------
// structure_file_in
// ----------------------------------------------------------------------------

/*!\brief A class for reading structured sequence files, e.g. Stockholm, Connect, Vienna, ViennaRNA bpp matrix ...
 * \ingroup structure_file
 * \tparam traits_type        An auxiliary type that defines certain member types and constants, must satisfy
 *                            seqan3::structure_file_input_traits_concept.
 * \tparam selected_field_ids A seqan3::fields type with the list and order of desired record entries; all fields
 *                            must be in seqan3::structure_file_in::field_ids.
 * \tparam valid_formats      A seqan3::type_list of the selectable formats (each must meet
 *                            seqan3::structure_file_input_format_concept).
 * \tparam stream_type        The type of the stream, must satisfy seqan3::istream_concept.
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
 * ```cpp
 * structure_file_in sf{"/tmp/my.dbn"}; // Vienna with RNA sequences assumed, regular std::ifstream taken as stream
 * ```
 *
 * Reading from an std::istringstream:
 *
 * ```cpp
 * std::string const input
 * {
 *     ">S.cerevisiae_tRNA-PHE M10740/1-73\n"
 *     "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA\n"
 *     "(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))). (-17.50)\n"
 *     "> example\n"
 *     "UUGGAGUACACAACCUGUACACUCUUUC\n"
 *     "..(((((..(((...)))..)))))... (-3.71)\n"
 * };
 *
 * std::istringstream iss(input);
 *
 * structure_file_in fin{std::move(iss), structure_file_format_vienna{}};
 * //              ^ no need to specify the template arguments
 * ```
 *
 * Note that this is not the same as writing `structure_file_in<>` (with angle brackets). In the latter case they are
 * explicitly set to their default values, in the former case
 * [automatic deduction](http://en.cppreference.com/w/cpp/language/class_template_argument_deduction) happens which
 * chooses different parameters depending on the constructor arguments. For opening from file, `structure_file_in<>`
 * would have also worked, but for opening from stream it would not have.
 *
 * In some cases, you do need to specify the arguments, e.g. if you want to read amino acids:
 *
 * ```cpp
 * structure_file_in<structure_file_default_traits_aa> fin{"/tmp/my.dbn"};
 * ```
 *
 * You can define your own traits type to further customise the types used by and returned by this class, see
 * seqan3::structure_file_default_traits_rna for more details. As mentioned above, specifying at least one
 * template parameter yourself means that you loose automatic deduction so if you want to read amino acids **and**
 * want to read from a string stream you need to give all types yourself:
 *
 * ```cpp
 *
 *  // ... input had amino acid sequences
 * std::istringstream iss(input);
 *
 * structure_file_in<structure_file_default_traits_aa,
 *                   fields<field::SEQ, field::ID, field::STRUCTURE>,
 *                   type_list<structure_file_format_vienna>,
 *                   std::istringstream> fin{std::move(iss), structure_file_format_vienna{}};
 * ```
 *
 * ### Reading record-wise
 *
 * You can iterate over this file record-wise:
 *
 * ```cpp
 * structure_file_in fin{"/tmp/my.dbn"};
 *
 * for (auto & rec : fin)
 * {
 *     std::cout << "ID: " << get<field::ID>(rec) << '\n';
 *     std::cout << "SEQ: " << (get<field::SEQ>(rec) | view::to_char) << '\n'; // sequence is converted to char on-the-fly
 *     std::cout << "STRUCTURE: " << (get<field::STRUCTURE>(rec) | view::to_char) << '\n';
 * }
 * ```
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
 * ```cpp
 * structure_file_in fin{"/tmp/my.dbn"};
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
 * structure_file_in fin{"/tmp/my.dbn"};
 *
 * for (auto & [ seq, id, structure ] : fin)
 * {
 *     std::cout << "ID: " << id << '\n';
 *     std::cout << "SEQ: " << (seq | view::to_char) << '\n'; // sequence is converted to char on-the-fly
 *     std::cout << "STRUCTURE: " << (structure | view::to_char) << '\n';
 * }
 * ```
 *
 * In this case you immediately get the three elements of the tuple: `seq` of \ref seq_type, `id` of
 * \ref id_type and `structure` of \ref structure_type. **But beware: with structured bindings you do need
 * to get the order of elements correctly!**
 *
 * ### Reading record-wise (custom fields)
 *
 * If you want to skip specific fields from the record you can pass a non-empty fields trait object to the
 * structure_file_in constructor to select the fields that should be read from the input. For example to choose a
 * combined field for SEQ and STRUCTURE (see above). Or to never actually read the STRUCTURE, if you don't need it.
 * The following snippets demonstrate the usage of such a fields trait object.
 *
 * ```cpp
 * structure_file_in fin{"/tmp/my.dbn", fields<field::ID, field::STRUCTURED_SEQ>{}};
 *
 * for (auto & [ id, structured_seq ] : fin) // note that the order is now different, "id" comes first, because it was specified first
 * {
 *     std::cout << "ID: " << id << '\n';
 *     // sequence and structure are part of the same vector, of type std::vector<structured_rna<rna5, wuss51>>
 *     std::cout << "SEQ: "  << (structured_seq | view::get<0> | view::to_char) << '\n'; // sequence string is extracted and converted to char on-the-fly
 *     std::cout << "STRUCTURE: " << (structured_seq | view::get<1> | view::to_char) << '\n'; // structure string is extracted and converted to char on-the-fly
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
 * structure_file_in fin{"/tmp/my.dbn"};
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
 *     concatenated_sequences<rna5_vector>         sequences;
 *     concatenated_sequences<std::string>         ids;
 *     concatenated_sequences<std::vector<wuss51>> structures;
 * };
 *
 * data_storage_t data_storage; // a global or globally used variable in your program
 *
 * // ... in your file reading function:
 *
 * structure_file_in fin{"/tmp/my.dbn"};
 *
 * data_storage.sequences = std::move(get<field::SEQ>(fin)); // we move the buffer directly into our storage
 * data_storage.ids = std::move(get<field::ID>(fin)); // we move the buffer directly into our storage
 * data_storage.structures = std::move(get<field::STRUCTURE>(fin)); // we move the buffer directly into our storage
 * ```
 *
 * Note that for this to make sense, your storage data types need to be identical to the corresponding column types
 * of the file. If you require different column types you can specify you own traits, see
 * seqan3::structure_file_input_traits_concept.
 *
 * ### Formats
 *
 * Currently, the only implemented format is seqan3::structure_file_format_vienna. More formats will follow soon.
 */
template<structure_file_input_traits_concept traits_type_ = structure_file_input_default_traits_rna,
         detail::fields_concept selected_field_ids_ = fields<field::SEQ, field::ID, field::STRUCTURE>,
         detail::type_list_of_structure_file_input_formats_concept valid_formats_
             = type_list<structure_file_format_vienna>,
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
    using iterator        = detail::in_file_iterator<structure_file_in>;
    //!\brief The const iterator type is void, because files are not const-iterable.
    using const_iterator  = void;
    //!\brief The type returned by end().
    using sentinel        = std::ranges::default_sentinel;
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
     * \param[in] _file_name Path to the file you wish to open.
     * \param[in] fields_tag A seqan3::fields tag. [optional]
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
                    if (std::ranges::equal(ext, extension))
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
     * \tparam file_format The format of the file in the stream, must satisfy
     * seqan3::structure_file_input_format_concept.
     * \param[in] _stream The stream to operate on (this must be std::move'd in!).
     * \param[in] format_tag The file format tag.
     * \param[in] fields_tag A seqan3::fields tag. [optional]
     */
    template<structure_file_input_format_concept file_format>
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
     * structure_file_in fin{"/tmp/my.dbn"};
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

    //!\brief Path of the file that the stream operates on.
    std::string file_name;

    //!\brief The stream we are reading from.
    stream_type stream;

    //!\brief File is at position 1 behind the last record.
    bool at_end{false};

    //!\brief Type of the format, an std::variant over the `valid_formats`.
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
        std::visit([&] (structure_file_input_format_concept & f)
        {
            // read new record
            if constexpr (selected_field_ids::contains(field::STRUCTURED_SEQ))
            {
                static_assert(!selected_field_ids::contains(field::STRUCTURE),
                              "You may not select field::STRUCTURED_SEQ and field::STRUCTURE at the same time.");
                static_assert(!selected_field_ids::contains(field::SEQ),
                              "You may not select field::STRUCTURED_SEQ and field::SEQ at the same time.");
                f.read(stream,
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
                f.read(stream,
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
 * \relates seqan3::structure_file_in
 * \{
 */
template <istream_concept<char> stream_type,
          structure_file_input_format_concept file_format,
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
/*!\brief std::tuple_size overload for column-like access.
 * [metafunction specialisation for seqan3::structure_file_in]
 */
template<seqan3::structure_file_input_traits_concept traits_type,
         seqan3::detail::fields_concept selected_field_ids,
         seqan3::detail::type_list_of_structure_file_input_formats_concept valid_formats,
         seqan3::istream_concept<char> stream_type>
struct tuple_size<seqan3::structure_file_in<traits_type, selected_field_ids, valid_formats, stream_type>>
{
    //!\brief The value equals the number of selected fields in the file.
    static constexpr size_t value = selected_field_ids::as_array.size();
};

/*!\brief std::tuple_element overload for column-like access.
 * [metafunction specialisation for seqan3::structure_file_in]
 */
template<size_t elem_no,
         seqan3::structure_file_input_traits_concept traits_type,
         seqan3::detail::fields_concept selected_field_ids,
         seqan3::detail::type_list_of_structure_file_input_formats_concept valid_formats,
         seqan3::istream_concept<char> stream_type>
struct tuple_element<elem_no, seqan3::structure_file_in<traits_type, selected_field_ids, valid_formats, stream_type>>
    : tuple_element<elem_no, typename seqan3::structure_file_in<traits_type,
                                                                selected_field_ids,
                                                                valid_formats,
                                                                stream_type>::file_as_tuple_type>
{};

} // namespace std

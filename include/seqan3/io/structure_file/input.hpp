// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::structure_file_input and corresponding traits classes.
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 */

#pragma once

#include <cassert>
#include <seqan3/std/filesystem>
#include <fstream>
#include <limits>
#include <optional>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include <seqan3/alphabet/adaptation/char.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/nucleotide/rna15.hpp>
#include <seqan3/alphabet/nucleotide/rna5.hpp>
#include <seqan3/alphabet/structure/dssp9.hpp>
#include <seqan3/alphabet/structure/structured_aa.hpp>
#include <seqan3/io/stream/concept.hpp>
#include <seqan3/io/exception.hpp>
#include <seqan3/io/detail/in_file_iterator.hpp>
#include <seqan3/io/detail/misc_input.hpp>
#include <seqan3/io/detail/record.hpp>
#include <seqan3/io/structure_file/input_format_concept.hpp>
#include <seqan3/io/structure_file/input_options.hpp>
#include <seqan3/io/structure_file/format_vienna.hpp>
#include <seqan3/io/structure_file/record.hpp>
#include <seqan3/utility/container/concept.hpp>
#include <seqan3/utility/type_list/traits.hpp>

namespace seqan3
{
// ----------------------------------------------------------------------------
// structure_file_input_traits
// ----------------------------------------------------------------------------

/*!\interface seqan3::structure_file_input_traits <>
 * \brief The requirements a traits_type for seqan3::structure_file_input must meet.
 * \ingroup io_structure_file
 */
/*!\name Requirements for seqan3::structure_file_input_traits
 * \brief You can expect these **member types** of all types that satisfy seqan3::structure_file_input_traits.
 * \memberof seqan3::structure_file_input_traits
 *
 * \details
 *
 * Note that the alphabet type of the seqan3::field::structured_seq is a combined alphabet, i.e. either
 * seqan3::structured_rna<sequence_alphabet, structure_alphabet> or
 * seqan3::structured_aa<sequence_alphabet, structure_alphabet> and the container type templates for
 * the field are those of seqan3::field::seq and seqan3::field::structure, respectively.
 *
 * \{
 */
/*!\typedef using seq_alphabet
 * \brief Alphabet of the characters for the seqan3::field::seq; must satisfy seqan3::alphabet.
 */
/*!\typedef using seq_legal_alphabet
 * \brief Intermediate alphabet for seqan3::field::seq; must satisfy seqan3::alphabet and be convertible to
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
 * \brief Type template of the seqan3::field::seq, a container template over `seq_alphabet`;
 * must satisfy seqan3::sequence_container.
 */
/*!\typedef using id_alphabet
 * \brief Alphabet of the characters for the seqan3::field::id; must satisfy seqan3::alphabet.
 */
/*!\typedef using id_container
 * \brief Type template of the seqan3::field::id, a container template over `id_alphabet`;
 * must satisfy seqan3::sequence_container.
 */
/*!\typedef using bpp_prob
 * \brief Data type for the base pair probabilities in seqan3::field::bpp; must satisfy std::is_floating_point.
 */
/*!\typedef using bpp_partner
 * \brief Data type for the partner index of an interaction in seqan3::field::bpp; must satisfy
 * std::numeric_limits::is_integer.
 */
/*!\typedef using bpp_queue
 * \brief A container template representing a set of interactions of type bpp_item,
 * which are (comparable) tuples of `bpp_prob` and `bpp_partner`;
 * must satisfy seqan3::container and must provide an std::emplace(bpp_prob, bpp_partner) function.
 */
/*!\typedef using bpp_container
 * \brief Type template of the seqan3::field::bpp, a container template over a set (bpp_queue) of interactions;
 * must satisfy seqan3::sequence_container.
 */
/*!\typedef using structure_alphabet
 * \brief Alphabet of the characters for the seqan3::field::structure; must satisfy seqan3::rna_structure_alphabet.
 */
/*!\typedef using structure_container
 * \brief Type template of the seqan3::field::structure, a container template over `structure_alphabet`;
 * must satisfy seqan3::sequence_container.
 */
/*!\typedef using energy_type
 * \brief Type template of the seqan3::field::energy; must be std::optional of a type satisfying std::is_floating_point.
 * \details
 * If the file record contains an energy, the value can be obtained through dereference or std::optional::value.
 * Otherwise, operator bool or std::optional::has_value return false.
 */
/*!\typedef using react_type
 * \brief Data type for the reactivity and reactivity error in seqan3::field::react and seqan3::field::react_err,
 * respectively; must satisfy std::is_floating_point.
 */
/*!\typedef using react_container
 * \brief Type template of the seqan3::field::react and seqan3::field::react_err, a container template over
 * `react_type`; must satisfy seqan3::sequence_container.
 */
/*!\typedef using comment_alphabet
 * \brief Alphabet of the characters for the seqan3::field::comment; must satisfy seqan3::alphabet.
 */
/*!\typedef using comment_container
 * \brief Type template of the seqan3::field::comment, a container template over `comment_alphabet`;
 * must satisfy seqan3::sequence_container.
 */
/*!\typedef using offset_type
 * \brief Type template of the seqan3::field::offset; must statisfy std::numeric_limits::is_integer.
 */
//!\}
//!\cond
template <typename t>
SEQAN3_CONCEPT structure_file_input_traits = requires(t v)
{
    // TODO(joergi-w) The expensive concept checks are currently omitted. Check again when compiler has improved.
    // sequence
    requires writable_alphabet<typename t::seq_alphabet>;
    requires writable_alphabet<typename t::seq_legal_alphabet>;
    requires explicitly_convertible_to<typename t::seq_legal_alphabet, typename t::seq_alphabet>;
    requires sequence_container<typename t::template seq_container<typename t::seq_alphabet>>;

    // id
    requires writable_alphabet<typename t::id_alphabet>;
    requires sequence_container<typename t::template id_container<typename t::id_alphabet>>;

    // bpp
    requires std::is_floating_point_v<typename t::bpp_prob>;
    requires std::numeric_limits<typename t::bpp_partner>::is_integer;

//    requires container // TODO check Associative container Concept when implemented
//        <typename t::template bpp_queue
//            <typename t::template bpp_item
//                <typename t::bpp_prob, typename t::bpp_partner>>>
//        && requires(typename t::template bpp_queue // TODO maybe implement also a version that allows emplace_back
//            <typename t::template bpp_item
//                <typename t::bpp_prob, typename t::bpp_partner>> value) { value.emplace(1.0, 1); };
//    requires sequence_container
//        <typename t::template bpp_container
//            <typename t::template bpp_queue
//                 <typename t::template bpp_item
//                      <typename t::bpp_prob, typename t::bpp_partner>>>>;

    // structure
    requires std::is_same_v<typename t::structure_alphabet, dssp9> // TODO(joergi-w) add aa_structure_concept
          || rna_structure_alphabet<typename t::structure_alphabet>;
    requires sequence_container<typename t::template structure_container<typename t::structure_alphabet>>;

    // structured sequence: tuple composites of seq and structure
    requires std::is_base_of_v<alphabet_tuple_base
        <typename t::template structured_seq_alphabet
            <typename t::seq_alphabet, typename t::structure_alphabet>,
             typename t::seq_alphabet, typename t::structure_alphabet>,
        typename t::template structured_seq_alphabet<typename t::seq_alphabet, typename t::structure_alphabet>>;
//    requires sequence_container
//        <typename t::template structured_seq_container
//            <typename t::template structured_seq_alphabet
//                <typename t::seq_alphabet, typename t::structure_alphabet>>>;

    // energy: std::optional of floating point number
    requires std::is_floating_point_v<typename t::energy_type::value_type>;

    // reactivity [error]
    requires std::is_floating_point_v<typename t::react_type>;
    requires sequence_container<typename t::template react_container<typename t::react_type>>;

    // comment
    requires writable_alphabet<typename t::comment_alphabet>;
    requires sequence_container<typename t::template comment_container<typename t::comment_alphabet>>;

    // offset
    requires std::numeric_limits<typename t::offset_type>::is_integer;
};
//!\endcond

// ----------------------------------------------------------------------------
// structure_file_input_default_traits
// ----------------------------------------------------------------------------

/*!\brief The default traits for seqan3::structure_file_input
 * \implements structure_file_input_traits
 * \ingroup io_structure_file
 *
 * \details
 *
 * If you wish to change a single or a few types from the default, just inherit from this class and
 * "overwrite" the respective type definitions.
 *
 * This example will make the file read into a smaller alphabet:
 *
 * \include test/snippet/io/structure_file/structure_file_input_mod_traits.cpp
 */
struct structure_file_input_default_traits_rna
{
    /*!\name Member types
     * \brief Definitions to satisfy seqan3::structure_file_input_traits.
     * \{
     */

    // sequence

    //!\brief The sequence alphabet is seqan3::rna5.
    using seq_alphabet                       = rna5;

    //!\brief The legal sequence alphabet for parsing is seqan3::rna15.
    using seq_legal_alphabet                 = rna15;

    //!\brief The type of an RNA sequence is std::vector.
    template <typename _seq_alphabet>
    using seq_container                      = std::vector<_seq_alphabet>;

    // id

    //!\brief The alphabet for an identifier string is char.
    using id_alphabet                        = char;

    //!\brief The string type for an identifier is std::basic_string.
    template <typename _id_alphabet>
    using id_container                       = std::basic_string<_id_alphabet>;

    // base pair probability structure

    //!\brief The type for a base pair probability is double.
    using bpp_prob                           = double;

    //!\brief The type for the partner position of a base pair probability is size_t.
    using bpp_partner                        = size_t;

    //!\brief The type of a base pair item is std::pair<double, size_t>.
    template <typename _bpp_prob, typename _bpp_partner>
    using bpp_item                           = std::pair<_bpp_prob, _bpp_partner>;

    //!\brief A queue of base pair items sorted by probability is realised with std::set.
    template <typename _bpp_item>
    using bpp_queue                          = std::set<_bpp_item>;

    //!\brief A string over all bases containing the respective interaction queues is represented as std::vector.
    template <typename _bpp_queue>
    using bpp_container                      = std::vector<_bpp_queue>;

    // fixed structure

    //!\brief The alphabet for a structure annotation is seqan3::phred42.
    using structure_alphabet                 = wuss51;

    //!\brief The string type for a structure annotation is std::vector.
    template <typename _structure_alphabet>
    using structure_container                = std::vector<_structure_alphabet>;

    // combined sequence and structure

    //!\brief The combined structured sequence alphabet is seqan3::structured_rna<seqan3::rna5, seqan3::wuss51>.
    template <typename _seq_alphabet, typename _structure_alphabet>
    using structured_seq_alphabet            = structured_rna<_seq_alphabet, _structure_alphabet>;

    //!\brief The type of a structured RNA sequence is std::vector.
    template <typename _structured_seq_alphabet>
    using structured_seq_container           = std::vector<_structured_seq_alphabet>;

    // energy

    //!\brief The type of the energy is std::optional<double>.
    using energy_type                        = std::optional<double>;

    // reactivity [error]

    //!\brief The type of the reactivity and reactivity error is double.
    using react_type                         = double;

    //!\brief The type of a string of reactivity values is std::vector.
    template <typename _react_type>
    using react_container                    = std::vector<_react_type>;

    // comment

    //!\brief The alphabet for a comment string is char.
    using comment_alphabet                   = char;

    //!\brief The string type for a comment is std::basic_string.
    template <typename _comment_alphabet>
    using comment_container                  = std::basic_string<_comment_alphabet>;

    // offset

    //!\brief The type of the offset is size_t.
    using offset_type                        = size_t;
    //!\}
};

//!\brief A traits type that specifies input as amino acids.
//!\ingroup io_structure_file
struct structure_file_input_default_traits_aa : structure_file_input_default_traits_rna
{
    /*!\name Member types
     * \brief Definitions to satisfy seqan3::structure_file_input_traits.
     * \{
     */

    //!\brief The sequence alphabet is seqan3::aa27.
    using seq_alphabet             = aa27;
    //!\brief The legal sequence alphabet for parsing is seqan3::aa27.
    using seq_legal_alphabet       = aa27;
    //!\brief The structure annotation alphabet is seqan3::dssp9.
    using structure_alphabet       = dssp9;
    //!\brief The combined structured sequence alphabet is seqan3::structured_aa<seqan3::aa27, seqan3::dssp9>.
    template <typename _seq_alphabet, typename _structure_alphabet>
    using structured_seq_alphabet  = structured_aa<_seq_alphabet, _structure_alphabet>;
    //!\}
};

// ----------------------------------------------------------------------------
// structure_file_input
// ----------------------------------------------------------------------------

/*!\brief A class for reading structured sequence files, e.g. Stockholm, Connect, Vienna, ViennaRNA bpp matrix ...
 * \ingroup io_structure_file
 * \tparam traits_type        An auxiliary type that defines certain member types and constants, must satisfy
 *                            seqan3::structure_file_input_traits.
 * \tparam selected_field_ids A seqan3::fields type with the list and order of desired record entries; all fields
 *                            must be in seqan3::structure_file_input::field_ids.
 * \tparam valid_formats      A seqan3::type_list of the selectable formats (each must meet
 *                            seqan3::structure_file_input_format).
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
 *   1. seqan3::field::seq (sequence)
 *   2. seqan3::field::id (identifier)
 *   3. seqan3::field::bpp (annotated sequence)
 *   4. seqan3::field::structure (secondary structure)
 *   5. seqan3::field::structured_seq (sequence and structure in one range)
 *   6. seqan3::field::energy (minimum free energy)
 *   7. seqan3::field::react (reactivity)
 *   8. seqan3::field::react_err (reactivity error)
 *   9. seqan3::field::comment (free text)
 *   10. seqan3::field::offset (index of first sequence character)
 *
 * The first three fields are retrieved by default (and in that order). The seqan3::field::structured_seq may be
 * selected to have sequence and structure directly stored in a more memory-efficient combined container.
 * If you select this field you must not select seqan3::field::seq or seqan3::field::structure.
 *
 * ### Construction and specialisation
 *
 * This class comes with two constructors, one for construction from a file name and one for construction from
 * an existing stream and a known format. The first one automatically picks the format based on the extension
 * of the file name. The second can be used if you have a non-file stream, like std::cin or std::istringstream,
 * that you want to read from and/or if you cannot use file-extension based detection, but know that your input
 * file has a certain format.
 *
 * In most cases the template parameters are deduced completely automatically, e.g. reading from a std::istringstream:
 * \include test/snippet/io/structure_file/structure_file_input_auto_temp_deduc.cpp
 *
 * Note that this is not the same as writing `structure_file_input<>` (with angle brackets). In the latter case they are
 * explicitly set to their default values, in the former case
 * [automatic deduction](https://en.cppreference.com/w/cpp/language/class_template_argument_deduction) happens which
 * chooses different parameters depending on the constructor arguments. For opening from file, `structure_file_input<>`
 * would have also worked, but for opening from stream it would not have.
 *
 * In some cases, you do need to specify the arguments, e.g. if you want to read amino acids:
 *
 * \include test/snippet/io/structure_file/structure_file_input_arg_spec.cpp
 *
 * You can define your own traits type to further customise the types used by and returned by this class, see
 * seqan3::structure_file_default_traits_rna for more details. As mentioned above, specifying at least one
 * template parameter yourself means that you loose automatic deduction so if you want to read amino acids **and**
 * want to read from a string stream you need to give all types yourself:
 *
 * \include test/snippet/io/structure_file/structure_file_input_trait_def.cpp
 *
 * ### Reading record-wise
 *
 * You can iterate over this file record-wise:
 *
 * \include test/snippet/io/structure_file/structure_file_input_record_iter.cpp
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
 * \include test/snippet/io/structure_file/structure_file_input_data_out.cpp
 *
 * ### Reading record-wise (decomposed records)
 *
 * Instead of using `get` on the record, you can also use
 * [structured bindings](https://en.cppreference.com/w/cpp/language/structured_binding)
 * to decompose the record into its elements:
 *
 * \include test/snippet/io/structure_file/structure_file_input_structured_bindings.cpp
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
 * \include test/snippet/io/structure_file/structure_file_input_skip_fields.cpp
 *
 * When reading a file, all fields not present in the file (but requested implicitly or via the `selected_field_ids`
 * parameter) are ignored.
 *
 * ### Views on files
 *
 * Since SeqAn files are ranges, you can also create views over files. A useful example is to filter the records
 * based on certain criteria, e.g. minimum length of the sequence field:
 *
 * \include test/snippet/io/structure_file/structure_file_input_filter_criteria.cpp
 *
 * ### End of file
 *
 * You can check whether a file is at end by comparing begin() and end() (if they are the same, the file is at end).
 *
 * ### Formats
 *
 * Currently, the only implemented format is seqan3::format_vienna. More formats will follow soon.
 */
template <structure_file_input_traits traits_type_ = structure_file_input_default_traits_rna,
          detail::fields_specialisation selected_field_ids_ = fields<field::seq, field::id, field::structure>,
          detail::type_list_of_structure_file_input_formats valid_formats_ = type_list<format_vienna>>
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
    //!\brief Character type of the stream(s).
    using stream_char_type   = char;
    //!\}

    /*!\brief The subset of seqan3::field IDs that are valid for this file; order corresponds to the types in
     * \ref field_types.
     */
    using field_ids = fields<field::seq,
                             field::id,
                             field::bpp,
                             field::structure,
                             field::structured_seq,
                             field::energy,
                             field::react,
                             field::react_err,
                             field::comment,
                             field::offset>;

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
                      return !(selected_field_ids::contains(field::structured_seq) &&
                               (selected_field_ids::contains(field::seq) ||
                               (selected_field_ids::contains(field::structure))));
                  }(), "You may not select field::structured_seq and either of field::seq and field::structure "
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
    using record_type    = structure_record<detail::select_types_with_ids_t<field_types, field_ids, selected_field_ids>,
                                            selected_field_ids>;
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
    using sentinel        = std::default_sentinel_t;
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
        primary_stream{new std::ifstream{}, stream_deleter_default}
    {
        primary_stream->rdbuf()->pubsetbuf(stream_buffer.data(), stream_buffer.size());
        static_cast<std::basic_ifstream<char> *>(primary_stream.get())->open(filename,
                                                                             std::ios_base::in | std::ios::binary);

        if (!primary_stream->good())
            throw file_open_error{"Could not open file " + filename.string() + " for reading."};

        // possibly add intermediate decompression stream
        secondary_stream = detail::make_secondary_istream(*primary_stream, filename);

        // initialise format handler
        detail::set_format(format, filename);
    }

    /*!\brief Construct from an existing stream and with specified format.
     * \tparam file_format The format of the file in the stream, must satisfy
     * seqan3::structure_file_input_format.
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
    template <input_stream stream_t, structure_file_input_format file_format>
    //!\cond
        requires std::same_as<typename std::remove_reference_t<stream_t>::char_type, char>
    //!\endcond
    structure_file_input(stream_t & stream,
                         file_format const & SEQAN3_DOXYGEN_ONLY(format_tag),
                         selected_field_ids const & SEQAN3_DOXYGEN_ONLY(fields_tag) = selected_field_ids{}) :
        primary_stream{&stream, stream_deleter_noop},
        format{detail::structure_file_input_format_exposer<file_format>{}}
    {
        static_assert(list_traits::contains<file_format, valid_formats>,
                      "You selected a format that is not in the valid_formats of this file.");

        // possibly add intermediate decompression stream
        secondary_stream = detail::make_secondary_istream(*primary_stream);
    }

    //!\overload
    template <input_stream stream_t, structure_file_input_format file_format>
    //!\cond
        requires std::same_as<typename std::remove_reference_t<stream_t>::char_type, char>
    //!\endcond
    structure_file_input(stream_t && stream,
                         file_format const & SEQAN3_DOXYGEN_ONLY(format_tag),
                         selected_field_ids const & SEQAN3_DOXYGEN_ONLY(fields_tag) = selected_field_ids{}) :
        primary_stream{new stream_t{std::move(stream)}, stream_deleter_default},
        format{detail::structure_file_input_format_exposer<file_format>{}}
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
     * \include test/snippet/io/structure_file/structure_file_input_ref_return.cpp
     *
     * It most situations using the iterator interface or a range-based for-loop are preferable to using front(),
     * because you can only move to the next record via the iterator.
     *
     * In any case, don't forget the reference! If you want to save the data from the record elsewhere, use move:
     *
     * \include test/snippet/io/structure_file/structure_file_input_move.cpp
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
    structure_file_input_options<typename traits_type::seq_legal_alphabet,
                                 selected_field_ids::contains(field::structured_seq)> options;

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

    //!\brief Type of the format, an std::variant over the `valid_formats`.
    using format_type = typename detail::variant_from_tags<valid_formats,
                                                           detail::structure_file_input_format_exposer>::type;
    //!\brief The actual std::variant holding a pointer to the detected/selected format.
    format_type format;
    //!\}

    //!\brief Tell the format to move to the next record and update the buffer.
    void read_next_record()
    {
        // clear the record
        record_buffer.clear();

        // store the current position in the position buffer
        position_buffer = secondary_stream->tellg();

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
            if constexpr (selected_field_ids::contains(field::structured_seq))
            {
                static_assert(!selected_field_ids::contains(field::structure),
                              "You may not select field::structured_seq and field::structure at the same time.");
                static_assert(!selected_field_ids::contains(field::seq),
                              "You may not select field::structured_seq and field::seq at the same time.");
                f.read_structure_record(*secondary_stream,
                                        options,
                                        detail::get_or_ignore<field::structured_seq>(record_buffer), // seq
                                        detail::get_or_ignore<field::id>(record_buffer),
                                        detail::get_or_ignore<field::bpp>(record_buffer),
                                        detail::get_or_ignore<field::structured_seq>(record_buffer), // structure
                                        detail::get_or_ignore<field::energy>(record_buffer),
                                        detail::get_or_ignore<field::react>(record_buffer),
                                        detail::get_or_ignore<field::react_err>(record_buffer),
                                        detail::get_or_ignore<field::comment>(record_buffer),
                                        detail::get_or_ignore<field::offset>(record_buffer));
            }
            else
            {
                f.read_structure_record(*secondary_stream,
                                        options,
                                        detail::get_or_ignore<field::seq>(record_buffer),
                                        detail::get_or_ignore<field::id>(record_buffer),
                                        detail::get_or_ignore<field::bpp>(record_buffer),
                                        detail::get_or_ignore<field::structure>(record_buffer),
                                        detail::get_or_ignore<field::energy>(record_buffer),
                                        detail::get_or_ignore<field::react>(record_buffer),
                                        detail::get_or_ignore<field::react_err>(record_buffer),
                                        detail::get_or_ignore<field::comment>(record_buffer),
                                        detail::get_or_ignore<field::offset>(record_buffer));
            }
        }, format);
    }

    //!\brief Befriend iterator so it can access the buffers.
    friend iterator;
};

/*!\name Type deduction guides
 * \relates seqan3::structure_file_input
 * \{
 */

//!\brief Deduction of the selected fields, the file format and the stream type.
template <input_stream stream_type,
          structure_file_input_format file_format,
          detail::fields_specialisation selected_field_ids>
structure_file_input(stream_type && stream, file_format const &, selected_field_ids const &)
    -> structure_file_input<typename structure_file_input<>::traits_type,       // actually use the default
                            selected_field_ids,
                            type_list<file_format>>;

//!\overload
template <input_stream stream_type,
          structure_file_input_format file_format,
          detail::fields_specialisation selected_field_ids>
structure_file_input(stream_type & stream, file_format const &, selected_field_ids const &)
    -> structure_file_input<typename structure_file_input<>::traits_type,       // actually use the default
                            selected_field_ids,
                            type_list<file_format>>;
//!\}

} // namespace seqan3

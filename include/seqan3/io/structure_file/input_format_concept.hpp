// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::structure_file_input_format.
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 */

#pragma once

#include <fstream>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include <seqan3/alphabet/nucleotide/rna5.hpp>
#include <seqan3/alphabet/structure/structured_rna.hpp>
#include <seqan3/alphabet/structure/wuss.hpp>
#include <seqan3/io/structure_file/input_options.hpp>
#include <seqan3/utility/type_list/type_list.hpp>

namespace seqan3::detail
{

/*!\brief Internal class used to expose the actual format interface to read structure records from the file.
 * \ingroup io_structure_file
 *
 * \tparam format_type The type of the format to be exposed.
 *
 * \details
 *
 * Exposes the protected member function `read_structure_record` from the given `format_type`, such that the file can
 * call the proper function for the selected format.
 */
template <typename format_type>
struct structure_file_input_format_exposer : public format_type
{
public:
    // Can't use `using format_type::read_structure_record` as it produces a hard failure in the format concept check
    // for types that do not model the format concept, i.e. don't offer the proper read_structure_record interface.
    //!\brief Forwards to the seqan3::structure_file_input_format::read_structure_record interface.
    template <typename... ts>
    void read_structure_record(ts &&... args)
    {
        format_type::read_structure_record(std::forward<ts>(args)...);
    }
};

} // namespace seqan3::detail

namespace seqan3
{

/*!\interface seqan3::structure_file_input_format <>
 * \brief The generic concept for structure file in formats.
 * \ingroup io_structure_file
 *
 * \details
 *
 * The details of this concept are only relevant to developers who wish to implement their own format.
 * The requirements for this concept are given as related functions and type traits.
 * Types that satisfy this concept are shown as "implementing this interface".
 *
 * \remark For a complete overview, take a look at \ref io_structure_file
 */
//!\cond
template <typename t>
concept structure_file_input_format = requires (detail::structure_file_input_format_exposer<t> & v,
                                                std::ifstream & f,
                                                structure_file_input_options<rna5, false> & options,
                                                rna5_vector & seq,
                                                std::string & id,
                                                std::vector<std::set<std::pair<double, size_t>>> & bpp,
                                                std::vector<wuss51> & structure,
                                                std::vector<structured_rna<rna5, wuss51>> & structured_seq,
                                                double energy,
                                                double react,
                                                double react_err,
                                                std::string & comment,
                                                size_t offset) {
    t::file_extensions;

    {
        v.read_structure_record(f, options, seq, id, bpp, structure, energy, react, react_err, comment, offset)
    } -> std::same_as<void>;

    {
        v.read_structure_record(f,
                                options,
                                seq,
                                id,
                                bpp,
                                std::ignore,
                                std::ignore,
                                std::ignore,
                                std::ignore,
                                std::ignore,
                                std::ignore)
    } -> std::same_as<void>;

    {
        v.read_structure_record(f,
                                options,
                                structured_seq,
                                id,
                                std::ignore,
                                structured_seq,
                                energy,
                                std::ignore,
                                std::ignore,
                                std::ignore,
                                std::ignore)
    } -> std::same_as<void>;

    // This one is required to be a valid call, but should always throw at run-time.
    {
        v.read_structure_record(f,
                                options,
                                std::ignore,
                                std::ignore,
                                std::ignore,
                                std::ignore,
                                std::ignore,
                                std::ignore,
                                std::ignore,
                                std::ignore,
                                std::ignore)
    } -> std::same_as<void>;
};
//!\endcond

// Workaround for https://github.com/doxygen/doxygen/issues/9379
#if SEQAN3_DOXYGEN_ONLY(1) 0
template <typename t>
class structure_file_input_format
{};
#endif

/*!\name Requirements for seqan3::structure_file_input_format
 * \brief You can expect these **members** on all types that implement seqan3::structure_file_input_format.
 * \memberof seqan3::structure_file_input_format
 * \{
 */
/*!\fn void read_structure_record(stream_type & stream,
 *                                structure_file_input_options<seq_legal_alph_type, structured_seq_combined> const & options,
 *                                seq_type & seq,
 *                                id_type & id,
 *                                bpp_type & bpp,
 *                                structure_type & structure,
 *                                energy_type & energy,
 *                                react_type & react,
 *                                react_type & react_err,
 *                                comment_type & comment,
 *                                offset_type & offset)
 * \brief Read from the specified stream and back-insert into the given field buffers.
 * \tparam stream_type      Input stream, must satisfy seqan3::input_stream_over with `char`.
 * \tparam seq_type         Type of the seqan3::field::seq input; must satisfy std::ranges::output_range
 * over a seqan3::alphabet.
 * \tparam id_type          Type of the seqan3::field::id input; must satisfy std::ranges::output_range
 * over a seqan3::alphabet.
 * \tparam bpp_type         Type of the seqan3::field::bpp input; must satisfy std::ranges::output_range
 * over a set of pair of types satisfying std::is_floating_point and std::numeric_limits::is_integer, respectively.
 * \tparam structure_type   Type of the seqan3::field::structure input; must satisfy std::ranges::output_range
 * over a seqan3::rna_structure_alphabet.
 * \tparam energy_type      Type of the seqan3::field::energy input; must satisfy std::is_floating_point.
 * \tparam react_type       Type of the seqan3::field::react and seqan3::field::react_err input;
 * must satisfy std::is_floating_point.
 * \tparam comment_type     Type of the seqan3::field::comment input; must satisfy std::ranges::output_range
 * over a seqan3::alphabet.
 * \tparam offset_type      Type of the seqan3::field::offset input; must satisfy std::numeric_limits::is_integer.
 * \param[in,out] stream    The input stream to read from.
 * \param[in]     options   File specific options passed to the format.
 * \param[out]    seq       The buffer for seqan3::field::seq input, i.e. the "sequence".
 * \param[out]    id        The buffer for seqan3::field::id input, e.g. the header line.
 * \param[out]    bpp       The buffer for seqan3::field::bpp input.
 * \param[out]    structure The buffer for seqan3::field::structure input.
 * \param[out]    energy    The buffer for seqan3::field::energy input.
 * \param[out]    react     The buffer for seqan3::field::react input.
 * \param[out]    react_err The buffer for seqan3::field::react_err input.
 * \param[out]    comment   The buffer for seqan3::field::comment input.
 * \param[out]    offset    The buffer for seqan3::field::offset input.
 *
 * \details
 *
 * ### Additional requirements
 *
 *   * The function must also accept std::ignore as parameter for any of the fields.
 *     [This is enforced by the concept checker!]
 *   * In this case the data read for that field shall be discarded by the format.
 *   * Instead of passing the fields seqan3::field::seq and seqan3::field::structure, you may also pass
 *     seqan3::field::structured_seq to both parameters. If you do, the std::ranges::range_value_t of the argument must be
 *     a specialisation of seqan3::structured_rna and the second template parameter to
 *     seqan3::structure_file_input_options must be set to true.
 */
/*!\var static inline std::vector<std::string> seqan3::structure_file_input_format::file_extensions
 * \brief The format type is required to provide a vector of all supported file extensions.
 */
//!\}

} // namespace seqan3

namespace seqan3::detail
{

/*!\brief Auxiliary value metafuncton that checks whether a type is a seqan3::type_list and all types meet
 * seqan3::structure_file_input_format [default is false].
 * \ingroup core
 * \see seqan3::type_list_specialisationOfstructure_file_input_formats
 */
template <typename t>
constexpr bool is_type_list_of_structure_file_input_formats_v = false;

/*!\brief Auxiliary value metafuncton that checks whether a type is a seqan3::type_list and all types meet
 * seqan3::structure_file_input_format [overload].
 * \ingroup core
 * \see seqan3::type_list_specialisationOfstructure_file_input_formats
 */
template <typename... ts>
constexpr bool is_type_list_of_structure_file_input_formats_v<type_list<ts...>> =
    (structure_file_input_format<ts> && ...);

/*!\brief Auxiliary concept that checks whether a type is a seqan3::type_list and all types meet
 * seqan3::structure_file_input_format.
 * \ingroup core
 * \see seqan3::is_type_list_of_structure_file_formats_v
 */
template <typename t>
concept type_list_of_structure_file_input_formats = is_type_list_of_structure_file_input_formats_v<t>;
} // namespace seqan3::detail

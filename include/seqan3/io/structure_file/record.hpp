// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::structure_record.
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 */

#pragma once

#include <seqan3/io/record.hpp>

namespace seqan3
{
/*!\brief The record type of seqan3::structure_file_input.
 * \ingroup io_structure_file
 * \implements seqan3::detail::record_like
 * \tparam field_types The types of the fields in this record as a seqan3::type_list.
 * \tparam field_ids A seqan3::fields type with seqan3::field IDs corresponding to field_types.
 *
 * \remark For a complete overview, take a look at \ref io_structure_file
 */
template <typename field_types, typename field_ids>
class structure_record : public record<field_types, field_ids>
{
    //!\brief The base class.
    using base_t = record<field_types, field_ids>;

    //!\brief The underlying std::tuple class.
    using tuple_base_t = typename base_t::base_type;

    //!\copydoc seqan3::record::field_constant
    template <field f>
    using field_constant = typename base_t::template field_constant<f>;

    using base_t::get_impl;

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    structure_record() = default;                                     //!< Defaulted.
    structure_record(structure_record const &) = default;             //!< Defaulted.
    structure_record & operator=(structure_record const &) = default; //!< Defaulted.
    structure_record(structure_record &&) = default;                  //!< Defaulted.
    structure_record & operator=(structure_record &&) = default;      //!< Defaulted.
    ~structure_record() = default;                                    //!< Defaulted.

    //!\brief Inherit std::tuple's and seqan3::record constructors.
    using base_t::base_t;
    //!\}

    //!\brief The identifier, usually a string.
    decltype(auto) id() &&
    {
        return get_impl(field_constant<seqan3::field::id>{}, static_cast<tuple_base_t &&>(*this));
    }
    //!\copydoc seqan3::structure_record::id
    decltype(auto) id() const &&
    {
        return get_impl(field_constant<seqan3::field::id>{}, static_cast<tuple_base_t const &&>(*this));
    }
    //!\copydoc seqan3::structure_record::id
    decltype(auto) id() &
    {
        return get_impl(field_constant<seqan3::field::id>{}, static_cast<tuple_base_t &>(*this));
    }
    //!\copydoc seqan3::structure_record::id
    decltype(auto) id() const &
    {
        return get_impl(field_constant<seqan3::field::id>{}, static_cast<tuple_base_t const &>(*this));
    }

    //!\brief The "sequence", usually a range of nucleotides or amino acids.
    decltype(auto) sequence() &&
    {
        return get_impl(field_constant<seqan3::field::seq>{}, static_cast<tuple_base_t &&>(*this));
    }
    //!\copydoc seqan3::structure_record::sequence
    decltype(auto) sequence() const &&
    {
        return get_impl(field_constant<seqan3::field::seq>{}, static_cast<tuple_base_t const &&>(*this));
    }
    //!\copydoc seqan3::structure_record::sequence
    decltype(auto) sequence() &
    {
        return get_impl(field_constant<seqan3::field::seq>{}, static_cast<tuple_base_t &>(*this));
    }
    //!\copydoc seqan3::structure_record::sequence
    decltype(auto) sequence() const &
    {
        return get_impl(field_constant<seqan3::field::seq>{}, static_cast<tuple_base_t const &>(*this));
    }

    //!\brief Fixed interactions, usually a string of structure alphabet characters.
    decltype(auto) sequence_structure() &&
    {
        return get_impl(field_constant<seqan3::field::structure>{}, static_cast<tuple_base_t &&>(*this));
    }
    //!\copydoc seqan3::structure_record::sequence_structure
    decltype(auto) sequence_structure() const &&
    {
        return get_impl(field_constant<seqan3::field::structure>{}, static_cast<tuple_base_t const &&>(*this));
    }
    //!\copydoc seqan3::structure_record::sequence_structure
    decltype(auto) sequence_structure() &
    {
        return get_impl(field_constant<seqan3::field::structure>{}, static_cast<tuple_base_t &>(*this));
    }
    //!\copydoc seqan3::structure_record::sequence_structure
    decltype(auto) sequence_structure() const &
    {
        return get_impl(field_constant<seqan3::field::structure>{}, static_cast<tuple_base_t const &>(*this));
    }

    //!\brief Energy of a folded sequence, represented by one float number.
    decltype(auto) energy() &&
    {
        return get_impl(field_constant<seqan3::field::energy>{}, static_cast<tuple_base_t &&>(*this));
    }
    //!\copydoc seqan3::structure_record::energy
    decltype(auto) energy() const &&
    {
        return get_impl(field_constant<seqan3::field::energy>{}, static_cast<tuple_base_t const &&>(*this));
    }
    //!\copydoc seqan3::structure_record::energy
    decltype(auto) energy() &
    {
        return get_impl(field_constant<seqan3::field::energy>{}, static_cast<tuple_base_t &>(*this));
    }
    //!\copydoc seqan3::structure_record::energy
    decltype(auto) energy() const &
    {
        return get_impl(field_constant<seqan3::field::energy>{}, static_cast<tuple_base_t const &>(*this));
    }

    //!\brief Base pair probability matrix of interactions, usually a matrix of float numbers.
    decltype(auto) base_pair_probability_matrix() &&
    {
        // this is computed
        return get_impl(field_constant<seqan3::field::bpp>{}, static_cast<tuple_base_t &&>(*this));
    }
    //!\copydoc seqan3::structure_record::base_pair_probability_matrix
    decltype(auto) base_pair_probability_matrix() const &&
    {
        return get_impl(field_constant<seqan3::field::bpp>{}, static_cast<tuple_base_t const &&>(*this));
    }
    //!\copydoc seqan3::structure_record::base_pair_probability_matrix
    decltype(auto) base_pair_probability_matrix() &
    {
        return get_impl(field_constant<seqan3::field::bpp>{}, static_cast<tuple_base_t &>(*this));
    }
    //!\copydoc seqan3::structure_record::base_pair_probability_matrix
    decltype(auto) base_pair_probability_matrix() const &
    {
        return get_impl(field_constant<seqan3::field::bpp>{}, static_cast<tuple_base_t const &>(*this));
    }

    // decltype(auto) reactivity(); // unused
    // decltype(auto) reactivity_errors(); // unused
    // decltype(auto) comment(); // unused
    // decltype(auto) base_qualities(); // unused
};
} // namespace seqan3

namespace std
{

/*!\brief Provides access to the number of elements in a tuple as a compile-time constant expression.
 * \implements seqan3::unary_type_trait
 * \relates seqan3::structure_record
 * \see std::tuple_size_v
 */
template <typename field_types, typename field_ids>
struct tuple_size<seqan3::structure_record<field_types, field_ids>> :
    tuple_size<typename seqan3::structure_record<field_types, field_ids>::base_type>
{};

/*!\brief Obtains the type of the specified element.
 * \implements seqan3::transformation_trait
 * \relates seqan3::structure_record
 * \see [std::tuple_element](https://en.cppreference.com/w/cpp/utility/tuple/tuple_element)
 */
template <size_t elem_no, typename field_types, typename field_ids>
struct tuple_element<elem_no, seqan3::structure_record<field_types, field_ids>> :
    tuple_element<elem_no, typename seqan3::structure_record<field_types, field_ids>::base_type>
{};

} // namespace std

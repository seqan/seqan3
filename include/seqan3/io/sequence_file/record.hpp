// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::sequence_record.
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 */

#pragma once

#include <seqan3/io/record.hpp>

namespace seqan3
{
/*!\brief The record type of seqan3::sequence_file_input.
 * \ingroup io_sequence_file
 * \implements seqan3::detail::record_like
 * \tparam field_types The types of the fields in this record as a seqan3::type_list.
 * \tparam field_ids A seqan3::fields type with seqan3::field IDs corresponding to field_types.
 *
 * \remark For a complete overview, take a look at \ref io_sequence_file
 */
template <typename field_types, typename field_ids>
class sequence_record : public record<field_types, field_ids>
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
    sequence_record() = default;                                    //!< Defaulted.
    sequence_record(sequence_record const &) = default;             //!< Defaulted.
    sequence_record & operator=(sequence_record const &) = default; //!< Defaulted.
    sequence_record(sequence_record &&) = default;                  //!< Defaulted.
    sequence_record & operator=(sequence_record &&) = default;      //!< Defaulted.
    ~sequence_record() = default;                                   //!< Defaulted.

    //!\brief Inherit std::tuple's and seqan3::record constructors.
    using base_t::base_t;
    //!\}

    //!\copybrief seqan3::field::id
    decltype(auto) id() &&
    {
        return get_impl(field_constant<seqan3::field::id>{}, static_cast<tuple_base_t &&>(*this));
    }
    //!\copydoc seqan3::sequence_record::id
    decltype(auto) id() const &&
    {
        return get_impl(field_constant<seqan3::field::id>{}, static_cast<tuple_base_t const &&>(*this));
    }
    //!\copydoc seqan3::sequence_record::id
    decltype(auto) id() &
    {
        return get_impl(field_constant<seqan3::field::id>{}, static_cast<tuple_base_t &>(*this));
    }
    //!\copydoc seqan3::sequence_record::id
    decltype(auto) id() const &
    {
        return get_impl(field_constant<seqan3::field::id>{}, static_cast<tuple_base_t const &>(*this));
    }

    //!\copybrief seqan3::field::seq
    decltype(auto) sequence() &&
    {
        return get_impl(field_constant<seqan3::field::seq>{}, static_cast<tuple_base_t &&>(*this));
    }
    //!\copydoc seqan3::sequence_record::sequence
    decltype(auto) sequence() const &&
    {
        return get_impl(field_constant<seqan3::field::seq>{}, static_cast<tuple_base_t const &&>(*this));
    }
    //!\copydoc seqan3::sequence_record::sequence
    decltype(auto) sequence() &
    {
        return get_impl(field_constant<seqan3::field::seq>{}, static_cast<tuple_base_t &>(*this));
    }
    //!\copydoc seqan3::sequence_record::sequence
    decltype(auto) sequence() const &
    {
        return get_impl(field_constant<seqan3::field::seq>{}, static_cast<tuple_base_t const &>(*this));
    }

    //!\copybrief seqan3::field::qual
    decltype(auto) base_qualities() &&
    {
        return get_impl(field_constant<seqan3::field::qual>{}, static_cast<tuple_base_t &&>(*this));
    }
    //!\copydoc seqan3::sequence_record::base_qualities
    decltype(auto) base_qualities() const &&
    {
        return get_impl(field_constant<seqan3::field::qual>{}, static_cast<tuple_base_t const &&>(*this));
    }
    //!\copydoc seqan3::sequence_record::base_qualities
    decltype(auto) base_qualities() &
    {
        return get_impl(field_constant<seqan3::field::qual>{}, static_cast<tuple_base_t &>(*this));
    }
    //!\copydoc seqan3::sequence_record::base_qualities
    decltype(auto) base_qualities() const &
    {
        return get_impl(field_constant<seqan3::field::qual>{}, static_cast<tuple_base_t const &>(*this));
    }
};
} // namespace seqan3

namespace std
{

/*!\brief Provides access to the number of elements in a tuple as a compile-time constant expression.
 * \implements seqan3::unary_type_trait
 * \relates seqan3::sequence_record
 * \see std::tuple_size_v
 */
template <typename field_types, typename field_ids>
struct tuple_size<seqan3::sequence_record<field_types, field_ids>> :
    tuple_size<typename seqan3::sequence_record<field_types, field_ids>::base_type>
{};

/*!\brief Obtains the type of the specified element.
 * \implements seqan3::transformation_trait
 * \relates seqan3::sequence_record
 * \see [std::tuple_element](https://en.cppreference.com/w/cpp/utility/tuple/tuple_element)
 */
template <size_t elem_no, typename field_types, typename field_ids>
struct tuple_element<elem_no, seqan3::sequence_record<field_types, field_ids>> :
    tuple_element<elem_no, typename seqan3::sequence_record<field_types, field_ids>::base_type>
{};

} // namespace std

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
 * \brief Provides the seqan3::record template and the seqan3::field enum.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <tuple>

#include <meta/meta.hpp>

#include <seqan3/core/type_list.hpp>
#include <seqan3/core/metafunction/template_inspection.hpp>

namespace seqan3
{

// ----------------------------------------------------------------------------
// enum field
// ----------------------------------------------------------------------------

/*!\brief An enumerator for the fields used in file formats.
 * \ingroup io
 *
 * Some of the fields are shared between formats.
 *
 * \todo Table which files take which fields
 */
enum class field
{
    SEQ,            //!< The "sequence", usually a range of nucleotides or amino acids.
    ID,             //!< The identifier, usually a string.
    QUAL,           //!< The qualities, usually in phred-score notation.
    SEQ_QUAL,       //!< Sequence and qualities combined in one range.
    // ...
    BPP,            //!< Base pair probability matrix of interactions, usually a matrix of float numbers.
    STRUCTURE,      //!< Fixed interactions, usually a string of structure alphabet characters.
    STRUCTURED_SEQ, //!< Sequence and fixed interactions combined in one range.
    ENERGY,         //!< Energy of a folded sequence, represented by one float number.
    REACT,          //!< Reactivity values of the sequence characters given in a vector of float numbers.
    REACT_ERR,      //!< Reactivity error values given in a vector corresponding to REACT.
    COMMENT,        //!< Comment field of arbitrary content, usually a string.
    OFFSET,         //!< Sequence start position, unsigned value.
    ALIGN,          //!< The alignment
    // ...
    USER_DEFINED_0, //!< Identifier for user defined file formats and specialisations.
    USER_DEFINED_1, //!< Identifier for user defined file formats and specialisations.
    USER_DEFINED_2, //!< Identifier for user defined file formats and specialisations.
    USER_DEFINED_3, //!< Identifier for user defined file formats and specialisations.
    USER_DEFINED_4, //!< Identifier for user defined file formats and specialisations.
    USER_DEFINED_5, //!< Identifier for user defined file formats and specialisations.
    USER_DEFINED_6, //!< Identifier for user defined file formats and specialisations.
    USER_DEFINED_7, //!< Identifier for user defined file formats and specialisations.
    USER_DEFINED_8, //!< Identifier for user defined file formats and specialisations.
    USER_DEFINED_9, //!< Identifier for user defined file formats and specialisations.
};

// ----------------------------------------------------------------------------
// fields
// ----------------------------------------------------------------------------

/*!\brief       A class template that holds a choice of seqan3::field.
 * \ingroup io
 * \tparam fs   The fields you wish to be present in the seqan3::record returned by your file.
 * \see         seqan3::record
 *
 * This class acts as a compile time list of seqan3::field elements. It is used in specialising file classes
 * to determine the elements in a seqan3::record.
 *
 * ### Example
 *
 * ```cpp
 * // specify custom field combination/order to file:
 * sequence_file_in fin{"/tmp/my.fasta", fields<field::ID, field::SEQ>{}};
 *
 * auto record = fin.front(); // get current record, in this case the first
 *
 * // record is tuple-like type, but allows access via field identifiers:
 * auto & id = get<field::ID>(record);
 * auto & seq = get<field::SEQ>(record);
 * ```
 *
 */
template <field ...fs>
struct fields
{
    //!\privatesection
    //!\brief The template parameters stored in an array for easy access.
    static constexpr std::array<field, sizeof...(fs)> as_array{fs...};

    //!\brief Special value that indicates that index_of() failed.
    static constexpr size_t npos = std::numeric_limits<size_t>::max();

    //!\brief Retrieve the position of field in the parameter pack.
    static constexpr size_t index_of(field f)
    {
        for (size_t i = 0; i < sizeof...(fs); ++i)
            if (as_array[i] == f)
                return i;
        return npos;
    }

    //!\brief Whether a field is contained in the parameter pack.
    static constexpr bool contains(field f)
    {
        return index_of(f) != npos;
    }

    static_assert([] () constexpr
                  {
                      for (size_t i = 0; i < as_array.size(); ++i)
                          for (size_t j = i + 1; j < as_array.size(); ++j)
                              if (as_array[i] == as_array[j])
                                  return false;

                      return true;
                  } (), "You may not include a field twice into fields<>.");
};

// ----------------------------------------------------------------------------
// record
// ----------------------------------------------------------------------------

/*!\brief The class template that file records are based on; behaves like an std::tuple.
 * \ingroup io
 * \implements seqan3::tuple_like_concept
 * \tparam field_types The types of the fields in this record as a seqan3::type_list.
 * \tparam field_ids   A seqan3::fields type with seqan3::field IDs corresponding to field_types.
 *
 * This class template behaves just like an std::tuple, with the exception that it provides an additional
 * get-interface that takes a seqan3::field identifier. The traditional get interfaces (via index and
 * via type) are also supported, but discouraged, because accessing via seqan3::field is unambiguous and
 * better readable.
 *
 * ### Example
 *
 * For input files this template is specialised automatically and provided by the file via its `record_type` member.
 * For output files you my define it locally and pass instances of this to the output file's `push_back()`.
 *
 * This is how it works:
 *
 * ```cpp
 *
 * using types        = type_list<dna4_vector, std::string, std::vector<phred42>>;
 * using types_as_ids =    fields<field::SEQ,  field::ID,   field::QUAL>;
 *
 * using record_type  = record<types, types_as_ids>;
 * // record_type now mimics std::tuple<std::string, dna4_vector, std::vector<phred42>>, the order also depends on selected_ids
 *
 * record_type my_record;
 * get<0>(my_record) = "the most important sequence in the database";   // access via index
 * get<field::SEQ>(my_record) = "ACGT"_dna4;                            // access via seqan3::field
 * get<std::string>(my_record) = "the least important sequence in the database";   // access via type
 * ```
 */
template <typename field_types, typename field_ids>
struct record : detail::transfer_template_args_onto_t<field_types, std::tuple>
{
public:
    //!\brief A specialisation of std::tuple.
    using base_type = detail::transfer_template_args_onto_t<field_types, std::tuple>;

    /*!\name Constructors, destructor and assignment
     * \brief Rule of five explicitly defaulted.
     * \{
     */
    record() = default;
    record(record const &) = default;
    record & operator=(record const &) = default;
    record(record &&) = default;
    record & operator=(record &&) = default;

    //!\brief Inherit std::tuple's constructors.
    using base_type::base_type;
    //!\}

    static_assert(field_types::size() == field_ids::as_array.size(),
                  "You must give as many IDs as types to seqan3::record.");

    //!\brief (Re-)Initialise all tuple elements with `{}`.
    void clear()
    {
        clear_impl(*this, std::make_integer_sequence<size_t, field_types::size()>{});
    }
private:
    //!\brief Auxiliary function for clear().
    template <size_t ...Is>
    constexpr void clear_impl(base_type & tup, std::integer_sequence<size_t, Is...> const &)
    {
        ((std::get<Is>(tup) = {}), ...);
    }
};

} // namespace seqan3

namespace std
{

/*!\brief Value metafunction specialisation for seqan3::record; returns number of elements in record.
 * \see seqan3::record
 */
template <typename field_types, typename field_ids>
struct tuple_size<seqan3::record<field_types, field_ids>>
{
    //!\brief The value member. Delegates to same value on base_type.
    static constexpr size_t value = tuple_size_v<typename seqan3::record<field_types, field_ids>::base_type>;
};

/*!\brief Value metafunction specialisation for seqan3::record; returns the type of an element in the record.
 * \see seqan3::record
 */
template <size_t elem_no, typename field_types, typename field_ids>
struct tuple_element<elem_no, seqan3::record<field_types, field_ids>>
{
    //!\brief The member type. Delegates to same type on base_type.
    using type = std::tuple_element_t<elem_no, typename seqan3::record<field_types, field_ids>::base_type>;
};

} // namespace std

namespace seqan3
{

/*!\name Free function get() interface for seqan3::record based on seqan3::field.
 * \brief This is the tuple interface via seqan3::field, e.g. `get<field::SEQ>(tuple)`.
 * \relates seqan3::record
 * \{
 */
template <field f, typename field_types, typename field_ids>
auto & get(record<field_types, field_ids> & r)
{
    static_assert(field_ids::contains(f), "The record does not contain the field you wish to retrieve.");
    return std::get<field_ids::index_of(f)>(r);
}

template <field f, typename field_types, typename field_ids>
auto const & get(record<field_types, field_ids> const & r)
{
    static_assert(field_ids::contains(f), "The record does not contain the field you wish to retrieve.");
    return std::get<field_ids::index_of(f)>(r);
}

template <field f, typename field_types, typename field_ids>
auto && get(record<field_types, field_ids> && r)
{
    static_assert(field_ids::contains(f), "The record does not contain the field you wish to retrieve.");
    return std::get<field_ids::index_of(f)>(std::move(r));
}

template <field f, typename field_types, typename field_ids>
auto const && get(record<field_types, field_ids> const && r)
{
    static_assert(field_ids::contains(f), "The record does not contain the field you wish to retrieve.");
    return std::get<field_ids::index_of(f)>(std::move(r));
}
//!\}

} // namespace seqan3

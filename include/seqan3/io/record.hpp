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

#include <seqan3/core/platform.hpp>

namespace seqan3
{

/*!\brief An enumerator for the fields used in file formats.
 * \ingroup io
 *
 * Some of the fields are shared between formats.
 */
enum class field
{
    SEQ,        //!< The "sequence", usually a range of nucleotides or amino acids.
    ID,         //!< The identifier, usually a string.
    QUAL,       //!< The qualities, usually in phred-score notation.
    SEQ_QUAL,   //!< Sequence and qualities combined in one range.
};

/*!\brief A class template that holds a choice of seqan3::field.
 * \ingroup io
 * \tparam fs The fields you wish to be present in the seqan3::record returned by your file.
 *
 * This class acts as a compile time list of seqan3::field elements. It used in specialising file classes
 * to determine the elements in a seqan3::records.
 *
 * TODO example
 */
template <field ...fs>
struct fields
{
//!\privatesection
    static constexpr std::array<field, sizeof...(fs)> as_array{fs...};
    //TODO static_assert that no types included twice


    static constexpr size_t index_of(field f)
    {
        for (size_t i = 0; i < sizeof...(fs); ++i)
            if (as_array[i] == f)
                return i;
        return std::numeric_limits<size_t>::max();
    }

    static constexpr bool contains(field f)
    {
        return index_of(f) != std::numeric_limits<size_t>::max();
    }
};

} // namespace seqan3

namespace seqan3::detail
{

/*!\brief A type metafunction that generates the tuple type for internal use in seqan3::record.
 * \ingroup io
 * \tparam field_types          The types of the fields available to the record in a meta::list.
 * \tparam field_types_as_ids   A seqan3::fields type with seqan3::field IDs corresponding to field_types.
 * \tparam selected_field_ids   A seqan3::fields type with the subset (and order) of the fields selected.
 * \tparam tup_t                The tuple type being aggregated.
 * \tparam field_no             The field we are currently processing.
 *
 * This metafunction "iterates" over the selected_field_ids and retrieves the corresponding typenames from
 * field_types via their identifer in field_types_as_ids. It recursively builds up tup_t which it returns
 * once the end of selected_field_ids is reached.
 */
template <typename field_types,
          typename field_types_as_ids,
          typename selected_field_ids,
          typename tup_t = std::tuple<>,
          size_t field_no = 0>
struct generate_record_tuple_type
{
    //!\brief The return type.
    using type = tup_t;
};

//!\cond
template <typename field_types,
          typename field_types_as_ids,
          typename selected_field_ids,
          typename tup_t,
          size_t field_no>
    requires field_no < selected_field_ids::as_array.size()
struct generate_record_tuple_type<field_types, field_types_as_ids, selected_field_ids, tup_t, field_no>
{
    using new_tup_t = decltype(
        std::tuple_cat(tup_t{},
                       std::tuple<meta::at_c<field_types,
                                             field_types_as_ids::index_of(selected_field_ids::as_array[field_no])>>{}));

    using type = typename generate_record_tuple_type<field_types,
                                                     field_types_as_ids,
                                                     selected_field_ids,
                                                     new_tup_t,
                                                     field_no + 1>::type; // iterator over the selected_fields
};
//!\endcond

} // namespace seqan3::detail


// template <size_t pos = 0, typename ...elemement_types>
// constexpr size_t get_position_in_tuple(field const target, std::tuple<elemement_types...> const & tup)
// {
//     if constexpr (pos == std::tuple_size_v<std::tuple<elemement_types...>>)
//         return pos;
//
//     if (std::get<pos>(tup) == target)
//         return pos;
//     else
//         return get_position_in_tuple<pos + 1>(target, tup);
// }

namespace seqan3
{

/*!\brief The class template that file record are based on; behaves like an std::tuple.
 * \ingroup io
 * \tparam field_types          The types of the fields available to this record in a meta::list.
 * \tparam field_types_as_ids   A seqan3::fields type with seqan3::field IDs corresponding to field_types.
 * \tparam selected_field_ids   A seqan3::fields type with the subset (and order) of the fields selected.
 *
 * This class behaves just like an std::tuple, with the exception that it provides an additional
 * get<>() interface that can be adressed not by absolute position or actual type of the element,
 * but also via a seqan3::field identifier.
 *
 * TODO example
 */
template <typename field_types,
          typename field_types_as_ids,
          typename selected_field_ids>
struct record
{
    //!\privatesection
//     using field_types_list = meta::list<field_types...>;
//     using selected_field_ids = selected_field_ids_;

    static_assert(selected_field_ids::as_array.size() <= field_types_as_ids::as_array.size()); //TODO add message

    //!\brief The type of the actual std::tuple. Is derived via detail::generate_record_tuple_type.
    using tuple_type = typename detail::generate_record_tuple_type<field_types, field_types_as_ids, selected_field_ids>::type;

    //!\brief An std::tuple as data member.
    tuple_type tuple;
};

/*!\brief A type metafunction analogue to std::tuple_element but with seqan3::field based selection.
 * \tparam f The seqan3::field you wish to query in the tuple-like type.
 * \tparam t The tuple-like type.
 *
 * This is the base template that is not defined. An overload currently only exists for seqan3::record.
 */
template <field f, typename t>
struct tuple_element;

/*!\brief A type metafunction shortcut for seqan3::tuple_element.
 * \tparam f The seqan3::field you wish to query in the tuple-like type.
 * \tparam t The tuple-like type.
 * \relates seqan3::tuple_element
 */
template <field f, typename t>
using tuple_element_t = typename tuple_element<f, t>::type;

/*!\brief A type metafunction analogue to std::tuple_element but with seqan3::field based selection.
 * [specialisation for seqan3::record]
 * \tparam f                    The seqan3::field you wish to query in your record.
 * \tparam field_types          The template parameters of seqan3::record.
 * \tparam field_types_as_ids   The template parameters of seqan3::record.
 * \tparam selected_field_ids   The template parameters of seqan3::record.
 * \relates seqan3::record
 */
template <field f, typename field_types, typename field_types_as_ids, typename selected_field_ids>
struct tuple_element<f, record<field_types, field_types_as_ids, selected_field_ids>>
{
    static_assert(selected_field_ids::contains(f), "The record does not contain the field you wish to retrieve.");
    using type = std::tuple_element_t<selected_field_ids::index_of(f),
                                     typename seqan3::record<field_types,
                                                             field_types_as_ids,
                                                             selected_field_ids>::tuple_type>;
};

}

//TODO does our stuff in namespace std already appear in doxygen?
namespace std
{

/*!\brief Value metafunction specialisation for seqan3::record; returns number of elements in tuple.
 * \see seqan3::record
 */
template <typename field_types, typename field_types_as_ids, typename selected_field_ids>
struct tuple_size<seqan3::record<field_types, field_types_as_ids, selected_field_ids>>
{
    static constexpr size_t value = tuple_size_v<typename seqan3::record<field_types, field_types_as_ids, selected_field_ids>::tuple_type>;
};

template <size_t elem_no, typename field_types, typename field_types_as_ids, typename selected_field_ids>
struct tuple_element<elem_no, seqan3::record<field_types, field_types_as_ids, selected_field_ids>>
{
    using type = tuple_element_t<elem_no, typename seqan3::record<field_types, field_types_as_ids, selected_field_ids>::tuple_type>;
};

} // namespace std

namespace seqan3
{

/*!\name Free function get() interface for seqan3::record based on index.
 * \brief This is the traditional tuple interface via position, e.g. `get<0>(tuple)`.
 * \relates seqan3::record
 * \{
 */
template <size_t i, typename field_types, typename field_types_as_ids, typename selected_field_ids>
auto & get(record<field_types, field_types_as_ids, selected_field_ids> & r)
{
    return std::get<i>(r.tuple);
}

template <size_t i, typename field_types, typename field_types_as_ids, typename selected_field_ids>
auto const & get(record<field_types, field_types_as_ids, selected_field_ids> const & r)
{
    return std::get<i>(r.tuple);
}

template <size_t i, typename field_types, typename field_types_as_ids, typename selected_field_ids>
auto && get(record<field_types, field_types_as_ids, selected_field_ids> && r)
{
    return std::get<i>(std::move(r.tuple));
}

template <size_t i, typename field_types, typename field_types_as_ids, typename selected_field_ids>
auto const && get(record<field_types, field_types_as_ids, selected_field_ids> const && r)
{
    return std::get<i>(std::move(r.tuple));
}
//!\}

/* get by type */
//TODO

/*!\name Free function get() interface for seqan3::record based on seqan3::field.
 * \brief This is the tuple interface via seqan3::field, e.g. `get<field::SEQ>(tuple)`.
 * \relates seqan3::record
 * \{
 */
template <field f, typename field_types, typename field_types_as_ids, typename selected_field_ids>
auto & get(record<field_types, field_types_as_ids, selected_field_ids> & r)
{
    static_assert(selected_field_ids::contains(f), "The record does not contain the field you wish to retrieve.");
    return std::get<selected_field_ids::index_of(f)>(r.tuple);
}

template <field f, typename field_types, typename field_types_as_ids, typename selected_field_ids>
auto const & get(record<field_types, field_types_as_ids, selected_field_ids> const & r)
{
    static_assert(selected_field_ids::contains(f), "The record does not contain the field you wish to retrieve.");
    return std::get<selected_field_ids::index_of(f)>(r.tuple);
}

template <field f, typename field_types, typename field_types_as_ids, typename selected_field_ids>
auto && get(record<field_types, field_types_as_ids, selected_field_ids> && r)
{
    static_assert(selected_field_ids::contains(f), "The record does not contain the field you wish to retrieve.");
    return std::get<selected_field_ids::index_of(f)>(std::move(r.tuple));
}

template <field f, typename field_types, typename field_types_as_ids, typename selected_field_ids>
auto const && get(record<field_types, field_types_as_ids, selected_field_ids> const && r)
{
    static_assert(selected_field_ids::contains(f), "The record does not contain the field you wish to retrieve.");
    return std::get<selected_field_ids::index_of(f)>(std::move(r.tuple));
}
//!\}

}

// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the seqan3::record template and the seqan3::field enum.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <tuple>

#include <meta/meta.hpp>

#include <seqan3/core/type_list.hpp>
#include <seqan3/core/type_traits/template_inspection.hpp>

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
 * The following table shows the usage of fields in the respective files
 * (Note that each valid format for a file must handle all of its fields):
 *
 * | Field          | Sequence IO | Alignment IO | Structure IO |
 * | -------------- | ----------- | ------------ | ------------ |
 * | SEQ            |      ✅      |      ✅       |       ✅      |
 * | ID             |      ✅      |      ✅       |       ✅      |
 * | QUAL           |      ✅      |      ✅       |       ✅      |
 * | SEQ_QUAL       |      ✅      |      ✅       |       ✅      |
 * | OFFSET         |             |      ✅       |       ✅      |
 * | BPP            |             |              |       ✅      |
 * | STRUCTURE      |             |              |       ✅      |
 * | STRUCTURED_SEQ |             |              |       ✅      |
 * | ENERGY         |             |              |       ✅      |
 * | REACT          |             |              |       ✅      |
 * | REACT_ERR      |             |              |       ✅      |
 * | COMMENT        |             |              |       ✅      |
 * | ALIGNMENT      |             |      ✅       |              |
 * | REF_ID         |             |      ✅       |              |
 * | REF_SEQ        |             |      ✅       |              |
 * | REF_OFFSET     |             |      ✅       |              |
 * | HEADER_PTR     |             |      ✅       |              |
 * | FLAG           |             |      ✅       |              |
 * | MATE           |             |      ✅       |              |
 * | MAPQ           |             |      ✅       |              |
 * | TAGS           |             |      ✅       |              |
 * | BIT_SCORE      |             |      ✅       |              |
 * | EVALUE         |             |      ✅       |              |
 */
enum class field
{
    // Fields used in multiple contexts ........................................
    SEQ,            //!< The "sequence", usually a range of nucleotides or amino acids.
    ID,             //!< The identifier, usually a string.
    QUAL,           //!< The qualities, usually in phred-score notation.
    SEQ_QUAL,       //!< Sequence and qualities combined in one range.
    OFFSET,         //!< Sequence (SEQ) relative start position (0-based), unsigned value.

    // Fields unique to structure io ...........................................
    BPP,            //!< Base pair probability matrix of interactions, usually a matrix of float numbers.
    STRUCTURE,      //!< Fixed interactions, usually a string of structure alphabet characters.
    STRUCTURED_SEQ, //!< Sequence and fixed interactions combined in one range.
    ENERGY,         //!< Energy of a folded sequence, represented by one float number.
    REACT,          //!< Reactivity values of the sequence characters given in a vector of float numbers.
    REACT_ERR,      //!< Reactivity error values given in a vector corresponding to REACT.
    COMMENT,        //!< Comment field of arbitrary content, usually a string.

    // Fields unique to alignment io ...........................................
    ALIGNMENT,      //!< The (pairwise) alignment stored in an seqan3::alignment object.
    REF_ID,         //!< The identifier of the (reference) sequence that SEQ was aligned to.
    REF_SEQ,        //!< The (reference) "sequence" information, usually a range of nucleotides or amino acids.
    REF_OFFSET,     //!< Sequence (REF_SEQ) relative start position (0-based), unsigned value.
    HEADER_PTR,     //!< A pointer to the seqan3::alignment_file_header object storing header information.

    FLAG,           //!< The alignment flag (bit information), `uint16_t` value.
    MATE,           //!< The mate pair information given as a std::tuple of reference name, offset and template length.
    MAPQ,           //!< The mapping quality of the SEQ alignment, usually a ohred-scaled score.
    TAGS,           //!< The optional tags in the SAM format, stored in a dictionary.

    BIT_SCORE,      //!< The bit score (statistical significance indicator), unsigned value.
    EVALUE,         //!< The e-value (length normalized bit score), `double` value.

    // User defined field aliases .. ...........................................
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
 * \snippet test/snippet/io/record.cpp usage_1
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
 * \implements seqan3::TupleLike
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
 * \snippet test/snippet/io/record.cpp usage_2
 */
template <typename field_types, typename field_ids>
struct record : detail::transfer_template_args_onto_t<field_types, std::tuple>
{
public:
    //!\brief A specialisation of std::tuple.
    using base_type = detail::transfer_template_args_onto_t<field_types, std::tuple>;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    record() = default;                           //!< Defaulted
    record(record const &) = default;             //!< Defaulted
    record & operator=(record const &) = default; //!< Defaulted
    record(record &&) = default;                  //!< Defaulted
    record & operator=(record &&) = default;      //!< Defaulted
    ~record() = default;                          //!< Defaulted

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

/*!\brief Provides access to the number of elements in a tuple as a compile-time constant expression.
 * \implements seqan3::UnaryTypeTrait
 * \relates seqan3::record
 * \see std::tuple_size_v
 */
template <typename field_types, typename field_ids>
struct tuple_size<seqan3::record<field_types, field_ids>>
{
    //!\brief The value member. Delegates to same value on base_type.
    static constexpr size_t value = tuple_size_v<typename seqan3::record<field_types, field_ids>::base_type>;
};

/*!\brief Obtains the type of the specified element.
 * \implements seqan3::TransformationTrait
 * \relates seqan3::record
 * \see [std::tuple_element](https://en.cppreference.com/w/cpp/utility/tuple/tuple_element)
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

//!\brief Free function get() for seqan3::record based on seqan3::field.
template <field f, typename field_types, typename field_ids>
auto & get(record<field_types, field_ids> & r)
{
    static_assert(field_ids::contains(f), "The record does not contain the field you wish to retrieve.");
    return std::get<field_ids::index_of(f)>(r);
}

//!\overload
template <field f, typename field_types, typename field_ids>
auto const & get(record<field_types, field_ids> const & r)
{
    static_assert(field_ids::contains(f), "The record does not contain the field you wish to retrieve.");
    return std::get<field_ids::index_of(f)>(r);
}

//!\overload
template <field f, typename field_types, typename field_ids>
auto && get(record<field_types, field_ids> && r)
{
    static_assert(field_ids::contains(f), "The record does not contain the field you wish to retrieve.");
    return std::get<field_ids::index_of(f)>(std::move(r));
}

//!\overload
template <field f, typename field_types, typename field_ids>
auto const && get(record<field_types, field_ids> const && r)
{
    static_assert(field_ids::contains(f), "The record does not contain the field you wish to retrieve.");
    return std::get<field_ids::index_of(f)>(std::move(r));
}
//!\}

} // namespace seqan3

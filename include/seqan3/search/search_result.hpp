// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::search_result.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <seqan3/std/concepts>
#include <exception>

#include <seqan3/core/detail/debug_stream_type.hpp>
#include <seqan3/core/detail/empty_type.hpp>
#include <seqan3/core/type_traits/template_inspection.hpp>

namespace seqan3::detail
{
struct policy_result_builder;
}

namespace seqan3
{

/*!\brief A data holder for the search result.
 * \ingroup search
 * \tparam query_id_type The type of the query_id.
 * \tparam cursor_type The type of the cursor.
 * \tparam reference_id_type The type of the reference_id.
 * \tparam reference_begin_pos_type The type of the reference_begin_pos.
 *
 * \if DEV
 * \note If the information is not available the type of the respective data is of seqan3::detail::empty_type.
 * \endif
 */
template <typename query_id_type, typename cursor_type, typename reference_id_type, typename reference_begin_pos_type>
class search_result
{
private:
    //!\brief Stores the query_id of the search result.
    query_id_type query_id_{};
    //!\brief Stores the cursor of the search result.
    cursor_type cursor_{};
    //!\brief Stores the reference_id of the search result.
    reference_id_type reference_id_{};
    //!\brief Stores the reference_begin_pos of the search result.
    reference_begin_pos_type reference_begin_pos_{};

    /*!\brief Construct from a query id and an index cursor.
     * \tparam query_id_type The type of the query_id.
     * \tparam cursor_type The type of the cursor.
     * \param[in] id The query id of the search result.
     * \param[in] cursor The cursor of the search result.
     */
    search_result(query_id_type id, cursor_type cursor) :
        query_id_(id), cursor_(std::move(cursor))
    {}

    /*!\brief Construct from a query id, a reference id and a begin position in the reference.
     * \tparam query_id_type The type of the query_id.
     * \tparam reference_id_type The type of the reference id.
     * \tparam reference_begin_pos_type The type of the reference begin position.
     * \param[in] q_id The query id of the search result.
     * \param[in] ref_id The reference id of the search result.
     * \param[in] ref_pos The reference begin position of the search result.
     */
    search_result(query_id_type q_id, reference_id_type ref_id, reference_begin_pos_type ref_pos) :
        query_id_(q_id), reference_id_(std::move(ref_id)), reference_begin_pos_(std::move(ref_pos))
    {}

    //!\cond
    // Grant the policy access to private constructors.
    friend detail::policy_result_builder;
    // Currently, the query id is set within the search result range. This needs to be adapted.
    template <typename search_algorithm_t, typename query_range_t>
    friend class search_result_range;
    //!\endcond

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    search_result() = default; //!< Defaulted.
    search_result(search_result const &) = default; //!< Defaulted.
    search_result(search_result &&) = default; //!< Defaulted.
    search_result & operator=(search_result const &) = default; //!< Defaulted.
    search_result & operator=(search_result &&) = default; //!< Defaulted.
    ~search_result() = default; //!< Defaulted.
    //!\}

    /*!\name Accessors
     * \brief Functions to access elements of the search result.
     * \{
     */
    //!\brief Returns the id of the query which produced this search result.
    constexpr query_id_type query_id() const noexcept
    {
        return query_id_;
    }

    /*!\brief Returns the index cursor pointing to the suffix array range where the query was found.
     * \sa seqan3::fm_index_cursor
     * \sa seqan3::bi_fm_index_cursor
     */
    constexpr cursor_type cursor() const noexcept(!(std::same_as<cursor_type, detail::empty_type>))
    {
        if constexpr (std::same_as<cursor_type, detail::empty_type>)
        {
            throw std::logic_error{"You tried to access the index cursor but it was not selected in the output "
                                   "configuration of the search."};
        }
        return cursor_;
    }

    /*!\brief Returns the reference id where the query was found.
     *
     * \details
     * The reference id is an arithmetic value that corresponds to the index of the reference text in the index.
     * The order is determined on construction of the index.
     */
    constexpr reference_id_type reference_id() const noexcept(!(std::same_as<reference_id_type, detail::empty_type>))
    {
        if constexpr (std::same_as<reference_id_type, detail::empty_type>)
        {
            throw std::logic_error{"You tried to access the reference id but it was not selected in the output "
                                   "configuration of the search."};
        }
        return reference_id_;
    }

    //!\brief Returns the reference begin positions where the query was found in the reference text (at `reference id`).
    constexpr reference_begin_pos_type reference_begin_pos() const
        noexcept(!(std::same_as<reference_begin_pos_type, detail::empty_type>))
    {
        if constexpr (std::same_as<reference_begin_pos_type, detail::empty_type>)
        {
            throw std::logic_error{"You tried to access the reference begin position but it was not selected in the "
                                   "output configuration of the search."};
        }
        return reference_begin_pos_;
    }
    //!\}

    /*!\name Comparison
     * \{
     */
    //!\brief Returns whether `lhs` and `rhs` are the same.
    friend bool operator==(search_result const & lhs, search_result const & rhs) noexcept
    {
        bool equality = lhs.query_id_ == rhs.query_id_;
        if constexpr (!std::is_same_v<cursor_type, detail::empty_type>)
            equality &= lhs.cursor_ == rhs.cursor_;
        if constexpr (!std::is_same_v<reference_id_type, detail::empty_type>)
            equality &= lhs.reference_id_ == rhs.reference_id_;
        if constexpr (!std::is_same_v<reference_begin_pos_type, detail::empty_type>)
            equality &= lhs.reference_begin_pos_ == rhs.reference_begin_pos_;

        return equality;
    }

    //!\brief Returns whether `lhs` and `rhs` are not the same.
    friend bool operator!=(search_result const & lhs, search_result const & rhs) noexcept
    {
        return !(lhs == rhs);
    }
    //!\}
};

/*!\brief Print the seqan3::search_result to seqan3::debug_stream.
 * \tparam char_t The underlying character type of the seqan3::debug_stream_type.
 * \tparam search_result_t A specialization of seqan3::search_result.
 * \param ds The stream.
 * \param res The search result to print.
 * \relates seqan3::debug_stream_type
 */
template <typename char_t, typename search_result_t>
//!\cond
    requires detail::is_type_specialisation_of_v<remove_cvref_t<search_result_t>, search_result>
//!\endcond
inline debug_stream_type<char_t> & operator<<(debug_stream_type<char_t> & ds, search_result_t && res)
{
    ds << "<query_id:" << res.query_id();
    if constexpr (!std::same_as<decltype(res.cursor()), detail::empty_type>)
        ds << ", cursor present";
    if constexpr (!std::same_as<decltype(res.reference_id()), detail::empty_type>)
        ds << ", reference_id:" << res.reference_id();
    if constexpr (!std::same_as<decltype(res.reference_begin_pos()), detail::empty_type>)
        ds << ", reference_pos:" << res.reference_begin_pos();
    ds << ">";

    return ds;
}

} // namespace seqan3

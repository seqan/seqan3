// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::search_result.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/type_traits/template_inspection.hpp>

namespace seqan3::detail
{

template <typename iterator_t,
          typename sequence_id_t = std::nullopt_t *,
          typename sequence_pos_t = std::nullopt_t *>
struct search_result_value_type
{
    size_t query_id{};
    iterator_t iterator{};
    sequence_id_t sequence_id{};
    sequence_pos_t sequence_pos{};

    search_result_value_type() = default;                                                //!< Defaulted
    search_result_value_type(search_result_value_type const &) = default;                //!< Defaulted
    search_result_value_type(search_result_value_type &&) = default;                     //!< Defaulted
    search_result_value_type & operator=(search_result_value_type const &) = default;    //!< Defaulted
    search_result_value_type & operator=(search_result_value_type &&) = default;         //!< Defaulted
    ~search_result_value_type() = default;                                               //!< Defaulted

    // No position.
    search_result_value_type(size_t _id, iterator_t _it) :
        query_id(_id), iterator(_it) {}

    // Only one size_t that indicates sequence position. sequence id is then 0.
    search_result_value_type(size_t _id, iterator_t _it, sequence_pos_t _pos) :
    query_id(_id), iterator(_it), sequence_id(0), sequence_pos(_pos) {}

    // Search returns a pair of sequence if and position.
    search_result_value_type(size_t _id, iterator_t _it, std::pair<sequence_id_t, sequence_pos_t> _pair) :
        query_id(_id), iterator(_it), sequence_id(_pair.first), sequence_pos(_pair.second) {}


};

} // namespace seqan3::detail

namespace seqan3
{

template <typename search_result_traits>
    requires detail::is_type_specialisation_of_v<search_result_traits, detail::search_result_value_type>
class search_result
{
private:
    search_result_traits data;

    using iterator_t          = decltype(data.iterator);
    using sequence_id_t       = decltype(data.sequence_id);
    using sequence_pos_t      = decltype(data.sequence_pos);

public:
    search_result(search_result_traits value) : data(value) {};

    search_result() = default;                                     //!< Defaulted
    search_result(search_result const &) = default;                //!< Defaulted
    search_result(search_result &&) = default;                     //!< Defaulted
    search_result & operator=(search_result const &) = default;    //!< Defaulted
    search_result & operator=(search_result &&) = default;         //!< Defaulted
    ~search_result() = default;                                    //!< Defaulted

    constexpr size_t query_id() const noexcept
    {
        return data.query_id;
    }

    constexpr iterator_t iterator() const noexcept
    {
        return data.iterator;
    }

    constexpr sequence_id_t sequence_id() const noexcept
    {
        static_assert(!std::is_same_v<sequence_id_t, std::nullopt_t *>,
                      "Identifier is not available but should.");
        return data.sequence_id;
    }

    constexpr sequence_pos_t sequence_pos() const noexcept
    {
        static_assert(!std::is_same_v<sequence_pos_t, std::nullopt_t *>,
                      "Identifier is not available but should.");
        return data.sequence_pos;
    }

};

} // namespace seqan3

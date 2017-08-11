// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2017, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2017, Knut Reinert & MPI Molekulare Genetik
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

/*!\cond DEV
 * \file
 * \ingroup view
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::view::join lazy specialisation.
 * \endcond
 */

#pragma once

#include <sdsl/int_vector.hpp>
#include <sdsl/rank_support_v5.hpp>
#include <sdsl/select_support_mcl.hpp>
#include <sdsl/sd_vector.hpp>
#include <sdsl/rrr_vector.hpp>

#include <seqan3/range/detail/random_access_iterator.hpp>
#include <seqan3/range/concept.hpp>

namespace seqan3::detail
{

// --------------------------------------------------------------------------
// class view_join_ra (lazy version of view_join_ra)
// --------------------------------------------------------------------------

/*!\brief The type returned by seqan3::view::join_ra if seqan3::view_join_flags::LAZY is set.
 * \implements seqan3::view_concept
 * \implements seqan3::random_access_range_concept
 * \tparam irng_t The reference type of the view.
 * \tparam flags_with_lazy_set seqan3::view_join_flags where seqan3::view_join_flags::LAZY is set
 * (seqan3::view_join_flags::SPARSE may or may not be).
 * \ingroup view
 *
 * The lazy specialisation performs no operations on creation and instead builds the support
 * data structures on-the-fly. Altogether the over-head is currently not worth it, as
 * the impaired access time of accessing 10% of all items is already greater than
 * construction time of the non-lazy implementation plus the non-impaired access time.
 *
 * \attention
 * This specialisation does not fulfil seqan3::sized_range_concept and it is not const-iterable, i.e.
 * iterating over the data structure changes its state. This also means that **all access needs
 * to be synchronised**.
 */
template <typename irng_t, seqan3::view_join_flags flags_with_lazy_set>
//!\cond
    requires random_access_range_concept<irng_t> &&
             random_access_range_concept<ranges::range_reference_t<irng_t>> &&
             sized_range_concept<ranges::range_reference_t<irng_t>> &&
             ((flags_with_lazy_set & seqan3::view_join_flags::LAZY) == seqan3::view_join_flags::LAZY)
//!\endcond
class view_join_ra<irng_t, flags_with_lazy_set>
{
public:
    //!\brief Expose the template parameter.
    static constexpr seqan3::view_join_flags flags = flags_with_lazy_set;

    /*!\name Member types
     * \{
     */
    //!\brief The reference_type (equals the template parameter `ref_t`).
    using reference         = ranges::range_reference_t<ranges::range_reference_t<irng_t>>;
    //!\brief The const_reference type.
    using const_reference   = reference;
    //!\brief The value_type (which equals the reference_type with any references removed.
    using value_type        = std::remove_cv_t<std::remove_reference_t<reference>>;
    //!\brief The iterator type of this view (a random access iterator).
//     using iterator          = random_access_iterator<view_join_ra const>;
    using iterator          = view_join_ra_iterator<view_join_ra>;
    //!\brief The const iterator type of this view (same as iterator, because it's a view).
    using const_iterator    = iterator;
    //!\brief The sentinel type of this view.
    using sentinel          = view_join_ra_sentinel<view_join_ra>;
    //!\brief An unsigned integer type, usually std::size_t.
    using size_type         = ranges::size_type_t<std::decay_t<ranges::range_reference_t<irng_t>>>;
    //!\brief A signed integer type, usually std::ptrdiff_t.
    using difference_type   = ranges::difference_type_t<std::decay_t<ranges::range_reference_t<irng_t>>>;
    //!\}

private:

    //!\brief Helper variable.
    static constexpr bool is_sparse =
        (flags_with_lazy_set & seqan3::view_join_flags::SPARSE) == seqan3::view_join_flags::SPARSE;

    //!\brief A data structure that caches rank and select to speed up repeated queries.
    struct pos_hint_t
    {
        //!\brief Cached rank.
        size_type rank;
        //!\brief Cached select, i.e. left interval border of rank-th sequence.
        size_type select;
        //!\brief Cached select of rank+1, i.e. right interval border of rank-th sequence.
        size_type select_next;
    };

    //!\brief A local variable for caching the position hint (available since there is no const interface anyway).
    pos_hint_t pos_hint{std::numeric_limits<size_type>::max(),
                        std::numeric_limits<size_type>::max(),
                        std::numeric_limits<size_type>::max()};

    /*!\name Friends
     * \brief Befriend the iterator and the sentinel so that they can access pos_hint_t.
     * \{
     */
    friend iterator;
    friend sentinel;
    //!\}

    //!\brief Type of the bit-vector.
    using bit_vector_t      = sdsl::sd_vector<>;
    //!\brief Type of the rank support data structure.
    using rank_support_t    = sdsl::rank_support_sd<>;
    //!\brief Type of the select support data structure.
    using select_support_t  = sdsl::select_support_sd<>;

    //!\brief Aggregation of the data members.
    struct data_t
    {
        //!\brief The input range (of ranges).
        irng_t & irange;
        //!\brief How far the bit-vector and support have been instantiated.
        size_type progress{0};
        //!\brief Whether we have indexed irange completely (progress == total number of elements).
        bool at_end{false};

        //!\brief A bit vector with the end-positions of underlying views marked as 1. [used when non-sparse]
        bit_vector_t end_positions = bit_vector_t();
        //!\brief Rank support for the bit vector. [used when non-sparse]
        rank_support_t rank_support{};
        //!\brief Select support for the bit vector. [used when non-sparse]
        select_support_t select_support{};

        //!\brief Vector of delimiters. [used when sparse]
        std::vector<size_type> delimiter{0u};
    };

    //!\brief All actual data is implicitly shared between copies.
    std::shared_ptr<data_t> data;

    //!\brief Fill the bitvector and build rank+select support.
    void init(size_type const until, pos_hint_t & hint)
    {
        if ((data->at_end) || (until < data->progress))
            return;

        // EXPERIMENTAL
        if constexpr (is_sparse)
        {
            auto it = ranges::begin(data->irange) + data->delimiter.size() - 1;
            auto it_end = ranges::end(data->irange);
            for (; (it != it_end) && (data->delimiter.back() <= until * 2); ++it)
                data->delimiter.push_back(data->delimiter.back() + ranges::size(*it));

            if (it == it_end)
                data->at_end = true;

            data->progress = data->delimiter.back();

            // adapt hint
            hint.rank        = data->delimiter.size() - 2;
            hint.select      = data->delimiter[hint.rank];
            hint.select_next = data->delimiter[hint.rank + 1];

            return;
        }
        else
        {
            size_type target_size = until * 2; // double the "resizes" to get amortised run-time down

            size_type total_length = 0;
            size_type count = 0;
            hint.rank = std::numeric_limits<decltype(hint.rank)>::max();
            for (auto it = ranges::begin(data->irange);
                target_size >= total_length;
                ++count, total_length += ranges::size(*(it++)))
            {
                //TODO SEQAN_UNLIKELY
                if ((total_length > until) && (hint.rank == std::numeric_limits<decltype(hint.rank)>::max()))
                {
                    // adapt hint
                    hint.rank = count - 1;
                    hint.select = total_length - ranges::size(data->irange[hint.rank]);
                    hint.select_next = total_length;
                }
                //TODO SEQAN_UNLIKELY
                if (it == ranges::end(data->irange))
                {
                    data->at_end = true;
                    target_size = total_length - 1;
                    break;
                }

            }
            assert(total_length > target_size);

            if (hint.rank == std::numeric_limits<decltype(hint.rank)>::max())
            {
                // adapt hint
                hint.rank = count - 1;
                hint.select = total_length - ranges::size(data->irange[hint.rank]);
                hint.select_next = total_length;
            }

            sdsl::bit_vector bivi(total_length + 1, 0, 1);

            std::size_t prev_size = 0;
            for (auto && elem : data->irange)
            {
                bivi[prev_size += ranges::size(elem)] = 1;
                if (prev_size == total_length)
                    break;
            }

            if constexpr (std::is_same_v<bit_vector_t, sdsl::bit_vector>)
            {
                std::swap(data->end_positions, bivi);
            } else // sd_vector doesn't have resize :'(
            {
                data->end_positions = bit_vector_t(bivi);
            }

            data->rank_support = rank_support_t(&data->end_positions);
            data->select_support = select_support_t(&data->end_positions);

            data->progress = target_size + 1;
        }
    }
public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    view_join_ra() = default;
    constexpr view_join_ra(view_join_ra const & rhs) = default;
    constexpr view_join_ra(view_join_ra && rhs) = default;
    constexpr view_join_ra & operator=(view_join_ra const & rhs) = default;
    constexpr view_join_ra & operator=(view_join_ra && rhs) = default;
    ~view_join_ra() = default;

    /*!\brief Construct from another range.
     * \param irange The input range (of ranges).
     */
    view_join_ra(irng_t & irange)
        : data{new data_t{irange}}
    {}
    //!\}

    /*!\name Iterators
     * \{
     */
    /*!\brief Returns an iterator to the first element of the container.
     * \returns Iterator to the first element.
     *
     * If the container is empty, the returned iterator will be equal to end().
     *
     * \attention Note that there are no const versions of begin() and end(); and
     * the non-const versions are **not thread safe**.
     *
     * \par Complexity
     *
     * Constant.
     *
     * \par Exceptions
     *
     * No-throw guarantee.
     */
    iterator begin() noexcept;

    /*!\brief Returns an iterator to the element following the last element of the container.
     * \returns Iterator to the first element.
     *
     * This element acts as a placeholder; attempting to dereference it results in undefined behaviour.
     *
     * \attention Note that there are no const versions of begin() and end(); and
     * the non-const versions are **not thread safe**.
     *
     * \par Complexity
     *
     * Constant.
     *
     * \par Exceptions
     *
     * No-throw guarantee.
     */
    sentinel end() noexcept;
    //!\}

    /*!\name Element access
     * \{
     */
    /*!\brief Return the i-th element.
     * \param i The element to retrieve.
     *
     * Caches a position hint and uses that on subsequent accesses.
     *
     * Accessing an element behind the last causes undefined behaviour. In debug mode an assertion checks the size of
     * the container.
     * \par Complexity
     *
     * * if i is in the hinted region: \f$O(1)\f$
     * * else, if not sparse: **amortised** \f$O(\log(^n/_m))\f$
     * * else, if sparse: **amortised** \f$O(\log(m))\f$
     *
     * \par Exceptions
     *
     * Basic exception guarantee.
     */
    reference operator[](size_type const i)
    {
        return (*this)[{i, pos_hint}];
    }

    /*!\brief Return the i-th element, but use the position hint and try to be faster.
     * \param hinted A pair of the to-be-retrieved position and the hint.
     *
     * If the accessed element is inside the same sub-range, this skips the rank and select queries.
     *
     * Accessing an element behind the last causes undefined behaviour. In debug mode an assertion checks the size of
     * the container.
     * \par Complexity
     *
     * * if i is in the hinted region: \f$O(1)\f$
     * * else, if not sparse: **amortised** \f$O(\log(^n/_m))\f$
     * * else, if sparse: **amortised** \f$O(\log(m))\f$
     *
     * \par Exceptions
     *
     * Basic exception guarantee.
     */
    reference operator[](std::tuple<size_type, pos_hint_t &> const & hinted)
    {
        auto & [i, hint] = hinted;
        init(i, hint);
        assert(i <= data->progress);

        if (i == hint.select_next) // we just incremented to the next subrange after hint
        {
            hint.rank++;
            hint.select = hint.select_next;
            hint.select_next += ranges::size(data->irange[hint.rank]);
        }
        else if ((i > hint.select_next) || (i < hint.select)) // we are in any other subrange than hint
        {
            // EXPERIMENTAL
            if constexpr (is_sparse)
            {
                size_type l = 0;
                size_type m = 0;
                size_type r = data->delimiter.size() - 1;

                while (l <= r)
                {
                    m = (l + r) / 2;

                    if (data->delimiter[m] > i)
                        r = m - 1;
                    else if (data->delimiter[m+1] > i)
                        break;
                    else
                        l = m + 1;
                }
                hint.rank        = m;
                hint.select      = data->delimiter[m];
                hint.select_next = data->delimiter[m+1];
            } else
            {
                hint.rank       = ((i + 1) < ranges::size(data->end_positions))
                                        ? data->rank_support(i + 1)
                                        : ranges::size(data->irange) - 1;
                assert(hint.rank < ranges::size(data->irange));

                hint.select = hint.rank ? data->select_support(hint.rank) : 0;
                assert(i >= hint.select);

                hint.select_next = hint.select + ranges::size(data->irange[hint.rank]);
                assert(i <= hint.select_next);
            }
        } // else the hint was correct

        return data->irange[hint.rank][i - hint.select];
    }
    /*!\brief Return the first element.
     *
     * Calling front on an empty container is undefined. In debug mode an assertion checks the size of the container.
     * \par Complexity
     *
     * Constant.
     *
     * \par Exceptions
     *
     * Basic exception guarantee.
     */
    reference front()
    {
        return (*this)[0];
    }
    //!\}

    //!\brief Implicit conversion to container types.
    template <typename container_type>
    operator container_type()
        requires sequence_concept<std::decay_t<container_type>> &&
                 std::is_same_v<std::decay_t<value_type>, std::decay_t<ranges::range_value_type_t<container_type>>>
    {
        container_type ret;
        for (auto && c : (*this))
            ret.push_back(c);
        return ret;
    }

    /*!\cond DEV
     * \brief Return the size of the support data structures in bytes.
     */
    size_type size_in_bytes() const
    {
        if constexpr (is_sparse)
            return data->delimiter.size() * sizeof(size_type);
        else
            return sdsl::size_in_bytes(data->end_positions) +
                   sdsl::size_in_bytes(data->rank_support) +
                   sdsl::size_in_bytes(data->select_support);
    }
    //!\endcond
};

// --------------------------------------------------------------------------
// class view_join_ra_iterator
// --------------------------------------------------------------------------

/*!\brief A custom iterator for seqan3::detail::view_join_ra that caches position hints.
 * \tparam view_join_ra_lazy_t A specialisation of seqan3::detail::view_join_ra where seqan3::view_join_flags::LAZY is
 * set.
 * \ingroup view
 */
template <typename view_join_ra_lazy_t>
//!\cond
    requires ((view_join_ra_lazy_t::flags & seqan3::view_join_flags::LAZY) == seqan3::view_join_flags::LAZY)
//!\endcond
class view_join_ra_iterator<view_join_ra_lazy_t> :
    public random_access_iterator_base<view_join_ra_lazy_t, view_join_ra_iterator<view_join_ra_lazy_t>>
{
private:
    //!\brief Shortcut.
    using base = random_access_iterator_base<view_join_ra_lazy_t, view_join_ra_iterator>;
    //!\brief Import private member types from parent.
    using typename base::position_type;

    /*!\name Data members
     * \brief Make the parent's private data members visible.
     * \{
     */
    using base::host;
    using base::pos;
    //!\}

    //!\brief The data structure that caches the position.
    typename view_join_ra_lazy_t::pos_hint_t hint{std::numeric_limits<typename view_join_ra_lazy_t::size_type>::max(),
                                                  std::numeric_limits<typename view_join_ra_lazy_t::size_type>::max(),
                                                  std::numeric_limits<typename view_join_ra_lazy_t::size_type>::max()};

public:
    /*!\name Member types
     * \brief Make the parent's member types visible.
     * \{
     */
    using typename base::difference_type;
    using typename base::value_type;
    using typename base::reference;
    using typename base::const_reference;
    using typename base::pointer;
    using typename base::iterator_category;
    //!\}

    //!\brief Import the parent's constructors.
    using base::base;

    //!\brief Constructor that takes a position hint in addition to host and pos.
    constexpr view_join_ra_iterator(view_join_ra_lazy_t & _host,
                                    position_type _pos,
                                    typename view_join_ra_lazy_t::pos_hint_t const & _hint) :
        view_join_ra_iterator{_host, _pos}
    {
        hint = _hint;
    }

    //!\brief Befriend the parent.
    friend base;

    //!\brief Befriend the sentinel.
    friend view_join_ra_sentinel<view_join_ra_lazy_t>;

    /*!\name Arithmetic operators
     * \{
     */
    //!\brief Pre-increment, return updated iterator.
    constexpr view_join_ra_iterator & operator++() noexcept
    {
        ++pos;
        host->init(pos, hint);
        return *this;
    }

    //!\brief Post-increment, return previous iterator state.
    constexpr view_join_ra_iterator operator++(int) noexcept
    {
        view_join_ra_iterator cpy{*this};
        ++pos;
        host->init(pos, hint);
        return cpy;
    }

    //!\brief Forward this iterator.
    constexpr view_join_ra_iterator & operator+=(difference_type const skip) noexcept
    {
        pos += skip;
        host->init(pos, hint);
        return *this;
    }

    //!\brief Forward copy of this iterator.
    constexpr view_join_ra_iterator operator+(difference_type const skip) const noexcept
    {
        view_join_ra_iterator cpy{*this};
        return cpy += skip;
    }

    //!\brief Decrement iterator by skip.
    constexpr view_join_ra_iterator & operator-=(difference_type const skip) noexcept
    {
        return *this += -skip;
    }

    //!\brief Return decremented copy of this iterator.
    constexpr view_join_ra_iterator operator-(difference_type const skip) const noexcept
    {
        return *this + (-skip);
    }

    //!\brief Return offset between this and remote iterator's position.
    constexpr difference_type operator-(view_join_ra_iterator const rhs) const noexcept
    {
        return pos - rhs.pos;
    }
    //!\}

    /*!\name Reference/Dereference operators
     * \brief The const-versions use the hint inside the view, while the non-const versions save
     * an iterator specific hint.
     * \{
    */
    //!\brief Dereference operator making use of the size hint.
    constexpr reference operator*() const noexcept(noexcept((*host)[{pos, hint}]))
    {
        return (*host)[pos];
    }

    //!\brief Return pointer to this iterator.
    constexpr pointer operator->() const noexcept(noexcept((&host)[{pos, hint}]))
    {
        return &host[pos];
    }

    //!\brief Return underlying container value currently pointed at.
    constexpr reference operator[](position_type const n) const noexcept(noexcept((*host)[{pos + n, hint}]))
    {
        return (*host)[pos + n];
    }
    //!\brief Dereference operator making use of the size hint.
    constexpr reference operator*() noexcept(noexcept((*host)[{pos, hint}]))
    {
        return (*host)[{pos, hint}];
    }

    //!\brief Return pointer to this iterator.
    constexpr pointer operator->() noexcept(noexcept((&host)[{pos, hint}]))
    {
        return &host[{pos, hint}];
    }

    //!\brief Return underlying container value currently pointed at.
    constexpr reference operator[](position_type const n) noexcept(noexcept((*host)[{pos + n, hint}]))
    {
        return (*host)[{pos + n, hint}];
    }
    //!\}
};

// --------------------------------------------------------------------------
// class view_join_ra_sentinel
// --------------------------------------------------------------------------

/*!\brief A custom sentinel for seqan3::detail::view_join_ra where seqan3::view_join_flags::LAZY is set.
 * \tparam view_join_ra_nolazy_t A specialisation of seqan3::detail::view_join_ra where seqan3::view_join_flags::LAZY is
 * set.
 * \ingroup view
 * \implements seqan3::sentinel_concept
 */
template <typename view_join_ra_lazy_t>
//!\cond
    requires ((view_join_ra_lazy_t::flags & seqan3::view_join_flags::LAZY) == seqan3::view_join_flags::LAZY)
//!\endcond
class view_join_ra_sentinel<view_join_ra_lazy_t>
{
private:
    //!\brief Pointer to host.
    view_join_ra_lazy_t * host = nullptr;

public:
    /*!\name Constructors/Destructors
     * \{
     */
    //!\brief Default constructor.
    constexpr view_join_ra_sentinel() = default;
    //!\brief Copy constructor.
    constexpr view_join_ra_sentinel(view_join_ra_sentinel const &) = default;
    //!\brief Copy construction via assignment.
    constexpr view_join_ra_sentinel & operator=(view_join_ra_sentinel const &) = default;
    //!\brief Move constructor.
    constexpr view_join_ra_sentinel (view_join_ra_sentinel &&) = default;
    //!\brief Move assignment.
    constexpr view_join_ra_sentinel & operator=(view_join_ra_sentinel &&) = default;
    //!\brief Use default deconstructor.
    ~view_join_ra_sentinel() = default;

    //!\brief Construct by host.
    explicit constexpr view_join_ra_sentinel(view_join_ra_lazy_t & host) noexcept : host{&host} {}
    //!\}

    /*!\name Comparison operators
     * \{
     */
    constexpr bool operator==(view_join_ra_iterator<view_join_ra_lazy_t> const & rhs) const noexcept
    {
        assert(rhs.host == host);
        return (rhs.host->data->at_end && (rhs.pos == rhs.host->data->progress));
    }

    constexpr bool operator!=(view_join_ra_iterator<view_join_ra_lazy_t> const & rhs) const noexcept
    {
        return !(*this == rhs);
    }

    constexpr bool operator<(view_join_ra_iterator<view_join_ra_lazy_t> const & rhs) const noexcept
    {
        assert(rhs.host == host);
        return false;
    }

    constexpr bool operator>(view_join_ra_iterator<view_join_ra_lazy_t> const & rhs) const noexcept
    {
        return !(*this == rhs);
    }

    constexpr bool operator<=(view_join_ra_iterator<view_join_ra_lazy_t> const & rhs) const noexcept
    {
        return (*this == rhs);
    }

    constexpr bool operator>=(view_join_ra_iterator<view_join_ra_lazy_t> const & rhs) const noexcept
    {
        assert(rhs.host == host);
        return true;
    }
    //!\}
};

/*!\name Comparison operators
 * \brief Compares only the absolute position of two iterators.
 * \relates seqan3::detail::view_join_ra_sentinel<view_join_ra_lazy_t>
 * \{
 */
template <typename view_join_ra_t>
constexpr bool operator==(view_join_ra_iterator<view_join_ra_t> const & lhs,
                          view_join_ra_sentinel<view_join_ra_t> const & rhs)
    requires ((view_join_ra_t::flags & seqan3::view_join_flags::LAZY) == seqan3::view_join_flags::LAZY)
{
    return (rhs == lhs);
}

template <typename view_join_ra_t>
constexpr bool operator!=(view_join_ra_iterator<view_join_ra_t> const & lhs,
                          view_join_ra_sentinel<view_join_ra_t> const & rhs)
    requires ((view_join_ra_t::flags & seqan3::view_join_flags::LAZY) == seqan3::view_join_flags::LAZY)
{
    return (rhs != lhs);
}

template <typename view_join_ra_t>
constexpr bool operator<(view_join_ra_iterator<view_join_ra_t> const & lhs,
                         view_join_ra_sentinel<view_join_ra_t> const & rhs)
    requires ((view_join_ra_t::flags & seqan3::view_join_flags::LAZY) == seqan3::view_join_flags::LAZY)
{
    return (rhs >= lhs);
}

template <typename view_join_ra_t>
constexpr bool operator>(view_join_ra_iterator<view_join_ra_t> const & lhs,
                         view_join_ra_sentinel<view_join_ra_t> const & rhs)
    requires ((view_join_ra_t::flags & seqan3::view_join_flags::LAZY) == seqan3::view_join_flags::LAZY)
{
    return (rhs <= lhs);
}

template <typename view_join_ra_t>
constexpr bool operator<=(view_join_ra_iterator<view_join_ra_t> const & lhs,
                          view_join_ra_sentinel<view_join_ra_t> const & rhs)
    requires ((view_join_ra_t::flags & seqan3::view_join_flags::LAZY) == seqan3::view_join_flags::LAZY)
{
    return (rhs > lhs);
}

template <typename view_join_ra_t>
constexpr bool operator>=(view_join_ra_iterator<view_join_ra_t> const & lhs,
                          view_join_ra_sentinel<view_join_ra_t> const & rhs)
    requires ((view_join_ra_t::flags & seqan3::view_join_flags::LAZY) == seqan3::view_join_flags::LAZY)
{
    return (rhs < lhs);
}
//!\}

// --------------------------------------------------------------------------
// class view_join_ra member definitions (post iterator definitions)
// --------------------------------------------------------------------------

template <typename irng_t, seqan3::view_join_flags flags_with_lazy_set>
//!\cond
    requires random_access_range_concept<irng_t> &&
             random_access_range_concept<ranges::range_reference_t<irng_t>> &&
             sized_range_concept<ranges::range_reference_t<irng_t>> &&
             ((flags_with_lazy_set & seqan3::view_join_flags::LAZY) == seqan3::view_join_flags::LAZY)
//!\endcond
inline typename view_join_ra<irng_t, flags_with_lazy_set>::iterator
view_join_ra<irng_t, flags_with_lazy_set>::begin() noexcept
{
    return typename view_join_ra<irng_t, flags_with_lazy_set>::iterator{
        *this, 0, {0, 0, ranges::size(data->irange[0])}};
}

template <typename irng_t, seqan3::view_join_flags flags_with_lazy_set>
//!\cond
    requires random_access_range_concept<irng_t> &&
             random_access_range_concept<ranges::range_reference_t<irng_t>> &&
             sized_range_concept<ranges::range_reference_t<irng_t>> &&
             ((flags_with_lazy_set & seqan3::view_join_flags::LAZY) == seqan3::view_join_flags::LAZY)
//!\endcond
inline typename view_join_ra<irng_t, flags_with_lazy_set>::sentinel
view_join_ra<irng_t, flags_with_lazy_set>::end() noexcept
{
    return typename view_join_ra<irng_t, flags_with_lazy_set>::sentinel{*this};
}

} // namespace seqan3::detail

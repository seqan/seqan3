// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::alignment_range.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <range/v3/utility/iterator.hpp>

#include <seqan3/std/concepts>
#include <seqan3/std/ranges>

namespace seqan3
{

/*!\brief The alignment
 * \ingroup execution
 * \tparam range_buffer_t The buffer type for the alignment stream.
 *
 * \details
 *
 * Provides a stream-like range interface over the alignments instances that are computed in a
 * seqan3::detail::alignment_executor_two_way executor.
 */
template <typename range_buffer_t>
class alignment_range
{
    static_assert(!std::is_const_v<range_buffer_t>,
                  "Cannot create an alignment stream over a const buffer.");

    //!\brief The iterator of seqan3::detail::alignment_range.
    class iterator_type
    {
    public:
        //!\brief Type for distances between iterators.
        using difference_type = typename alignment_range::difference_type;
        //!\brief Value type of container elements.
        using value_type = typename alignment_range::value_type;
        //!\brief Use reference type defined by container.
        using reference = typename alignment_range::reference;
        //!\brief Use const reference type provided by container.
        using const_reference = std::add_const_t<reference>;
        //!\brief Pointer type is pointer of container element type.
        using pointer = std::add_pointer_t<value_type>;
        //!\brief Sets iterator category as input iterator.
        using iterator_category = std::input_iterator_tag;

        /*!\name Constructors, destructor and assignment
         * \{
         */
        iterator_type()                                  = default;
        iterator_type(iterator_type const &)             = default;
        iterator_type(iterator_type &&)                  = default;
        iterator_type & operator=(iterator_type const &) = default;
        iterator_type & operator=(iterator_type &&)      = default;
        ~iterator_type()                                 = default;

        //!\brief Construct from alignment stream.
        iterator_type(alignment_range & _stream) : stream_ptr(&_stream)
        {}
        //!}

        /*!\name Read
         * \{
         */
        reference operator*()
        {
            return stream_ptr->cached();
        }

        const_reference operator*() const
        {
            return stream_ptr->cached();
        }
        //!\}

        /*!\name Increment operators
         * \{
         */
        iterator_type & operator++(/*pre*/)
        {
            stream_ptr->next();
            return *this;
        }

        iterator_type operator++(int /*post*/)
        {
            auto tmp{*this};
            stream_ptr->next();
            return tmp;
        }
        //!\}

        /*!\name Comparison operators
         * \{
         */

        constexpr bool operator==(std::ranges::default_sentinel const &) const
        {
            return stream_ptr->eof();
        }

        friend constexpr bool operator==(std::ranges::default_sentinel const & lhs,
                                         iterator_type const & rhs)
        {
            return rhs == lhs;
        }

        constexpr bool operator!=(std::ranges::default_sentinel const & rhs) const
        {
            return !(*this == rhs);
        }

        friend constexpr bool operator!=(std::ranges::default_sentinel const & lhs,
                                         iterator_type const & rhs)
        {
            return rhs != lhs;
        }
        //!\}
    private:
        //!\brief Pointer to the underlying range.
        alignment_range * stream_ptr{};
    };

    // Befriend the iterator with this class.
    friend class iterator_type;

public:

    //!\brief The offset type.
    using difference_type = typename range_buffer_t::difference_type;
    //!\brief The alignment result type.
    using value_type      = typename range_buffer_t::value_type;
    //!\brief The reference type.
    using reference       = typename range_buffer_t::reference;
    //!\brief The iterator type.
    using iterator        = iterator_type;
    //!\brief The sentinel type.
    using sentinel        = std::ranges::default_sentinel;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    alignment_range()                                    = delete;
    alignment_range(alignment_range const &)             = default;
    alignment_range(alignment_range &&)                  = default;
    alignment_range & operator=(alignment_range const &) = default;
    alignment_range & operator=(alignment_range &&)      = default;
    ~alignment_range()                                   = default;

    explicit alignment_range(range_buffer_t & _range_buffer) :
        range_buffer{&_range_buffer, [](auto) { /*no-op*/ }}
    {}

    // Construct from resource.
    template <typename resource_t,
              typename kernel_t>
    //!\cond
        requires std::is_constructible_v<range_buffer_t, resource_t, kernel_t>
    //!\endcond
    alignment_range(resource_t _range_buffer_resource,
                    kernel_t _kernel) :
        range_buffer{std::make_shared<range_buffer_t>(std::forward<resource_t>(_range_buffer_resource), _kernel)}
    {}
    //!}

    /*!\name Iterators
     * \{
     */
    iterator begin()
    {
        next();
        return iterator{*this};
    }

    sentinel end() noexcept
    {
        return {};
    }
    //!\}

protected:

    //!\brief Receives the next alignment result from the executor buffer.
    void next()
    {
        assert(!eof());
        if (auto opt = range_buffer->bump(); opt.has_value())
            cache = &(*opt).get();
        else
            eof_flag = true;
    }

    //!\brief Returns the cached result.
    // TODO make iterable in case of multiple results per alignment.
    auto & cached() noexcept
    {
        assert(cache);
        return *cache;
    }

    //!\copydoc cached()
    auto const & cached() const noexcept
    {
        assert(cache);
        return *cache;
    }

    //!\brief Returns whether the executor buffer reached is end.
    constexpr bool eof() const noexcept
    {
        return eof_flag;
    }

private:
    //!\brief The underlying executor buffer.
    std::shared_ptr<range_buffer_t> range_buffer;
    //!\brief Stores last read element.
    value_type * cache{};
    //!\brief Indicates whether the stream has reached its end.
    bool eof_flag{false};
};

} // namespace seqan3

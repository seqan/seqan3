// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
//
// Copyright (c) 2006-2017, Knut Reinert, FU Berlin
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
// ==========================================================================

#pragma once

#include <iterator>
#include <type_traits>

#include <range/v3/utility/iterator.hpp>

#include <seqan3/core/concept/iterator.hpp>
#include <seqan3/range/container/concept.hpp>
#include <seqan3/range/concept.hpp>

#include <seqan3/io/detail/chunking.hpp>

namespace seqan3::detail
{

template <typename container_t>
class chunk_adaptor_iterator_base
{
public:
    using difference_type   = typename container_t::difference_type; // -> global meta-functions.
    using value_type        = typename container_t::value_type; // -> global meta-function.
    using reference         = typename container_t::reference;  // -> global meta-function.
    using const_reference   = typename container_t::const_reference;    // -> global meta-function
    using pointer           = std::add_pointer_t<value_type>;

protected:

    using iterator_t = std::conditional_t<std::is_const_v<container_t>,
                                          typename container_t::const_iterator,
                                          typename container_t::iterator>;

    chunk_adaptor_iterator_base() = default;

    // custom constructor.
    chunk_adaptor_iterator_base(container_t & c, bool _to_end = false) :
        chunk_b(begin(c)),
        chunk_c(begin(c)),
        chunk_e(end(c)),
        cont_ptr(&c)
    {
        if (_to_end)
            chunk_c = chunk_e;
    }

    chunk_adaptor_iterator_base(chunk_adaptor_iterator_base const & /*other*/) = default;
    chunk_adaptor_iterator_base(chunk_adaptor_iterator_base && /*other*/) = default;

    chunk_adaptor_iterator_base & operator=(chunk_adaptor_iterator_base const & /*other*/) = default;
    chunk_adaptor_iterator_base & operator=(chunk_adaptor_iterator_base && /*other*/) = default;

    ~chunk_adaptor_iterator_base() = default;

    iterator_t chunk_b{};
    iterator_t chunk_c{};
    iterator_t chunk_e{};

    container_t * cont_ptr = nullptr;

    // member functions.

    inline iterator_t chunk_current() const
    {
        return chunk_c;
    }

    inline iterator_t chunk_end() const
    {
        return chunk_e;
    }

    template <typename integral_t>
    inline void
    advance_chunk_impl(integral_t const offset)
    {
        using std::advance;
        advance(chunk_c, offset);
    }

    inline void
    trim_trailing_impl() noexcept
    {
        // no-op
    }
};

template <typename container_t>
    requires input_range_concept<container_t>
class input_chunk_adaptor_iterator : public chunk_adaptor_iterator_base<container_t>,
                             public chunk_decorator<input_chunk_adaptor_iterator<container_t>>
{
    friend class chunk_decorator<input_chunk_adaptor_iterator<container_t>>;

    using base_iter_t = chunk_adaptor_iterator_base<container_t>;
    using chunk_t     = chunk_decorator<input_chunk_adaptor_iterator>;

public:

    using iterator_category = std::input_iterator_tag;

    input_chunk_adaptor_iterator() = default;

    // custom iterator.
    input_chunk_adaptor_iterator(container_t & c, bool _to_end = false) :
        base_iter_t(c, _to_end)
    {}

    input_chunk_adaptor_iterator(input_chunk_adaptor_iterator const & /*other*/) = default;
    input_chunk_adaptor_iterator(input_chunk_adaptor_iterator && /*other*/) = default;

    input_chunk_adaptor_iterator & operator=(input_chunk_adaptor_iterator const & /*other*/) = default;
    input_chunk_adaptor_iterator & operator=(input_chunk_adaptor_iterator && /*other*/) = default;

    ~input_chunk_adaptor_iterator() = default;

    typename base_iter_t::const_reference operator*() const
    {
        return *this->chunk_c;
    }

    input_chunk_adaptor_iterator & operator++(/*pre*/)
    {
        ++this->chunk_c;
        return *this;
    }

    input_chunk_adaptor_iterator operator++(int /*post*/)
    {
        auto tmp{*this};
        ++(*this);
        return tmp;
    }

    inline bool operator==(input_chunk_adaptor_iterator const & rhs) const
    {
        return this->chunk_c == rhs.chunk_c;
    }

    inline bool operator!=(input_chunk_adaptor_iterator const & rhs) const
    {
        return !(*this == rhs) ;
    }

private:

    template <typename integral_t>
    inline void
    next_chunk_impl(integral_t const /*chunk_size*/)
    {}
};

template <typename container_t>
    requires random_access_sequence_concept<container_t>
class output_chunk_adaptor_iterator : public chunk_adaptor_iterator_base<container_t>,
                              public chunk_decorator<output_chunk_adaptor_iterator<container_t>>
{
    friend class chunk_decorator<output_chunk_adaptor_iterator<container_t>>;
    using base_iter_t = chunk_adaptor_iterator_base<container_t>;
    using chunk_t = chunk_decorator<output_chunk_adaptor_iterator>;

public:

    using iterator_category = std::output_iterator_tag;

    output_chunk_adaptor_iterator() = default;

    // custom iterator.
    explicit output_chunk_adaptor_iterator(container_t & c) :
        base_iter_t(c, true)
    {}

    output_chunk_adaptor_iterator(output_chunk_adaptor_iterator const & /*other*/) = default;
    output_chunk_adaptor_iterator(output_chunk_adaptor_iterator && /*other*/) = default;

    output_chunk_adaptor_iterator & operator=(output_chunk_adaptor_iterator const & /*other*/) = default;
    output_chunk_adaptor_iterator & operator=(output_chunk_adaptor_iterator && /*other*/) = default;

    ~output_chunk_adaptor_iterator() = default;

    // operator= implementation.
    template <typename value_t>
    output_chunk_adaptor_iterator & operator=(value_t const & val)
    {
        if (this->chunk_c == this->chunk_e)  // appends one value.
        {
            next_chunk_impl(1);
        }
        *this->chunk_c = val;
        ++this->chunk_c;
        return *this;
    }

    // operator implementation.
    output_chunk_adaptor_iterator & operator*()
    {
        return *this;
    }

    output_chunk_adaptor_iterator & operator++(/*pre*/)
    {
        return *this;
    }

    output_chunk_adaptor_iterator operator++(int /*post*/)
    {
        return *this;
    }

    template <typename rhs_iterator>
    inline bool operator==(rhs_iterator const & rhs)
        requires sentinel_concept<rhs_iterator, output_chunk_adaptor_iterator>
    {
        return this->chunk_c == rhs.chunk_c;
    }

    template <typename rhs_iterator>
    inline bool operator!=(rhs_iterator const & rhs)
        requires sentinel_concept<rhs_iterator, output_chunk_adaptor_iterator>
    {
        return !(*this == rhs) ;
    }

private:

    template <typename integral_t>
    inline void
    next_chunk_impl(integral_t const chunk_size)
    {
        if (this->chunk_c == this->chunk_e)
        {
            auto pos = this->chunk_c - this->chunk_b;
            // call resize member function of pointed-to container.
            this->cont_ptr->resize(std::size(*this->cont_ptr) + chunk_size);
            // update pointer if capacity was not enough.
            this->chunk_b = std::begin(*this->cont_ptr);
            this->chunk_c = this->chunk_b + pos;
            this->chunk_e = std::end(*this->cont_ptr);
        }
    }

    inline void
    trim_trailing_impl()
    {
        this->cont_ptr->resize(this->chunk_c - this->chunk_b);
        this->chunk_e = std::end(*this->cont_ptr);
    }
};

template <typename container_t>
inline auto
make_preferred_input_iterator_range(container_t & c)
    requires input_range_concept<std::decay_t<container_t>>
{
    return std::tuple{input_chunk_adaptor_iterator{c}, input_chunk_adaptor_iterator{c, true}};
}

template <typename container_t>
inline auto
make_preferred_output_iterator(container_t & c)
    requires random_access_sequence_concept<std::decay_t<container_t>>
{
    return output_chunk_adaptor_iterator{c};
}

} //namespace seqan3::detail

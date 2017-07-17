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

#include <iterator>
#include <type_traits>

#include <range/v3/utility/iterator.hpp>

#include <seqan3/range/container/concept.hpp>
#include <seqan3/core/meta/associated_types.hpp>
#include <seqan3/io/detail/chunking.hpp>

namespace seqan3::detail
{

template <typename container_t>
class chunk_iterator_container_adaptor
{
public:
    // TODO(rrahn): make to global meta_function
    using difference_type   = typename container_t::difference_type; // -> global meta-functions.
    using value_type        = typename container_t::value_type; // -> global meta-function.
    using reference         = typename container_t::reference;  // -> global meta-function.
    using const_reference   = typename container_t::const_reference;    // -> global meta-function
    using pointer           = std::add_pointer_t<value_type>;

protected:

    using iterator_t = iterator_type_t<container_t>;  // global meta-function.

    chunk_iterator_container_adaptor() = default;

    // custom constructor.
    chunk_iterator_container_adaptor(container_t & c, bool _to_end = false) :
        chunk_b(begin(c)),
        chunk_c(begin(c)),
        chunk_e(end(c)),
        cont_ptr(&c)
    {
        if (_to_end)
            chunk_c = chunk_e;
    }

    chunk_iterator_container_adaptor(chunk_iterator_container_adaptor const & /*other*/) = default;
    chunk_iterator_container_adaptor(chunk_iterator_container_adaptor && /*other*/) = default;

    chunk_iterator_container_adaptor & operator=(chunk_iterator_container_adaptor const & /*other*/) = default;
    chunk_iterator_container_adaptor & operator=(chunk_iterator_container_adaptor && /*other*/) = default;

    ~chunk_iterator_container_adaptor() = default;

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
};

template <typename container_t>
//   requires resizable_container_concept<container_t>
class chunk_input_iterator : public chunk_iterator_container_adaptor<container_t>,
                             public chunk_decorator<chunk_input_iterator<container_t>>
{
    friend class chunk_decorator<chunk_input_iterator<container_t>>;

    using base_iter_t = chunk_iterator_container_adaptor<container_t>;
    using chunk_t     = chunk_decorator<chunk_input_iterator>;

public:

    // using base_iter_t::value_type;
    // using base_iter_t::difference_type;
    // using base_iter_t::reference;
    // using base_iter_t::const_reference;
    // using base_iter_t::pointer;
    using iterator_category = std::input_iterator_tag;

    chunk_input_iterator() = default;

    // custom iterator.
    chunk_input_iterator(container_t & c, bool _to_end = false) :
        base_iter_t(c, _to_end)
    {}

    chunk_input_iterator(chunk_input_iterator const & /*other*/) = default;
    chunk_input_iterator(chunk_input_iterator && /*other*/) = default;

    chunk_input_iterator & operator=(chunk_input_iterator const & /*other*/) = default;
    chunk_input_iterator & operator=(chunk_input_iterator && /*other*/) = default;

    ~chunk_input_iterator() = default;

    typename base_iter_t::const_reference operator*() const
    {
        return *this->chunk_c;
    }

    chunk_input_iterator & operator++(/*pre*/)
    {
        ++this->chunk_c;
        return *this;
    }

    chunk_input_iterator operator++(int /*post*/)
    {
        auto tmp{*this};
        ++(*this);
        return tmp;
    }

    inline bool
    _equal(chunk_input_iterator const & other) const
    {
        return this->chunk_c == other.chunk_c;
    }

private:
    // strong exception guarantee
    template <typename integral_t>
    inline void
    next_chunk_impl(integral_t const /*chunk_size*/)
    {}
};

template <typename container_t>
inline bool
operator==(chunk_input_iterator<container_t> const & lhs,
           chunk_input_iterator<container_t> const & rhs)
{
    return lhs._equal(rhs);
}

template <typename container_t>
inline bool
operator!=(chunk_input_iterator<container_t> const & lhs,
           chunk_input_iterator<container_t> const & rhs)
{
    return !lhs._equal(rhs);
}

template <typename container_t>
//   requires resizable_container_concept<container_t>
class chunk_output_iterator : public chunk_iterator_container_adaptor<container_t>,
                              public chunk_decorator<chunk_output_iterator<container_t>>
{
    friend class chunk_decorator<chunk_output_iterator<container_t>>;
    using base_iter_t = chunk_iterator_container_adaptor<container_t>;
    using chunk_t = chunk_decorator<chunk_output_iterator>;

public:

    // using base_iter_t::value_type;
    // using base_iter_t::difference_type;
    // using base_iter_t::reference;
    // using base_iter_t::const_reference;
    // using base_iter_t::pointer;
    using iterator_category = std::output_iterator_tag;

    chunk_output_iterator() = default;

    // custom iterator.
    explicit chunk_output_iterator(container_t & c) :
        base_iter_t(c, true)
    {}

    chunk_output_iterator(chunk_output_iterator const & /*other*/) = default;
    chunk_output_iterator(chunk_output_iterator && /*other*/) = default;

    chunk_output_iterator & operator=(chunk_output_iterator const & /*other*/) = default;
    chunk_output_iterator & operator=(chunk_output_iterator && /*other*/) = default;

    ~chunk_output_iterator() = default;

    // operator= implementation.
    template <typename value_t>
//        requires std::is_convertible_to_concept<value_type, value_t>
    chunk_output_iterator & operator=(value_t const & val)
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
    chunk_output_iterator & operator*()
    {
        return *this;
    }

    chunk_output_iterator & operator++(/*pre*/)
    {
        return *this;
    }

    chunk_output_iterator operator++(int /*post*/)
    {
        return *this;
    }

private:
    // strong exception guarantee
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
};

// Returns input itertor pointing to begin of container.
template <typename container_t>
//    requires container_concept<container_t>
inline auto
input_iterator(container_t & c)
{
    return ranges::v3::make_iterator_range(chunk_input_iterator{c}, chunk_input_iterator{c, true});
}

// What would be the output semantic? -> Appending to the container.
// What if we have a source buffer *? -> Start at the beginning of that pointer.
// Returns back_insert iterator for containers supporting *.push_back().
template <typename container_t>
//    requires sequence_concept<container_t>
inline auto
output_iterator(container_t & c)
{
    return chunk_output_iterator{c};
}

} //namespace seqan3::detail

namespace std
{

}

#ifdef NDEBUG

#include <string>

namespace seqan3::detail
{

// static_assert(std::is_same_v<value_type_t<chunk_input_iterator<std::string>>, char>);

}  // namespace seqan3::detail

#endif

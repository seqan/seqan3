// -*- C++ -*-
//===------------------------------ span ---------------------------------===//
//
//                     The LLVM Compiler Infrastructure
//
// This file is dual licensed under the MIT and the University of Illinois Open
// Source Licenses. See LICENSE.TXT for details.
//
//===---------------------------------------------------------------------===//


/*!\file
// \brief Provides std::span from the C++20 standard library.
// \see https://en.cppreference.com/w/cpp/container/span
*/

//!\cond
#pragma once

#include <array>        // for array
#include <cstddef>      // for ptrdiff_t
#include <iterator>     // for iterators
#include <type_traits>  // for remove_cv, etc

#include <seqan3/core/platform.hpp>

#if __has_include(<span>)
#include <span>
#else

namespace std
{
inline constexpr ptrdiff_t dynamic_extent = -1;
template <typename span_tp, ptrdiff_t span_extent = dynamic_extent> class span;

template <class span_tp>
struct is_span_impl : public false_type {};

template <class span_tp, ptrdiff_t span_extent>
struct is_span_impl<span<span_tp, span_extent>> : public true_type {};

template <class span_tp>
struct is_span : public is_span_impl<remove_cv_t<span_tp>> {};

template <class span_tp>
struct is_std_array_impl : public false_type {};

template <class span_tp, size_t span_sz>
struct is_std_array_impl<array<span_tp, span_sz>> : public true_type {};

template <class span_tp>
struct is_std_array : public is_std_array_impl<remove_cv_t<span_tp>> {};

template <class span_tp, class ElementType, class = void>
struct is_span_compatible_container : public false_type {};

template <class span_tp, class ElementType>
struct is_span_compatible_container<span_tp, ElementType,
        void_t<
        // is not a specialization of span
            typename enable_if<!is_span<span_tp>::value, nullptr_t>::type,
        // is not a specialization of array
            typename enable_if<!is_std_array<span_tp>::value, nullptr_t>::type,
        // is_array_v<Container> is false,
            typename enable_if<!is_array_v<span_tp>, nullptr_t>::type,
        // data(cont) and size(cont) are well formed
            decltype(data(declval<span_tp>())),
            decltype(size(declval<span_tp>())),
        // remove_pointer_t<decltype(data(cont))>(*)[] is convertible to ElementType(*)[]
            typename enable_if<
                is_convertible_v<remove_pointer_t<decltype(data(declval<span_tp &>()))>(*)[],
                                 ElementType(*)[]>,
                nullptr_t>::type
        >>
    : public true_type {};

template <typename span_tp, ptrdiff_t span_extent>
class span {
public:
//  constants and types
    using element_type           = span_tp;
    using value_type             = remove_cv_t<span_tp>;
    using index_type             = size_t;
    using difference_type        = ptrdiff_t;
    using pointer                = span_tp *;
    using const_pointer          = const span_tp *; // not in standard
    using reference              = span_tp &;
    using const_reference        = const span_tp &; // not in standard
    using iterator               = value_type *;
    using const_iterator         = value_type const *;
    using reverse_iterator       = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;

    static constexpr index_type extent = span_extent;
    static_assert(span_extent >= 0, "Can't have a span with an extent < 0");

// [span.cons], span constructors, copy, assignment, and destructor
    constexpr span() noexcept : data{nullptr}
    { static_assert(span_extent == 0, "Can't default construct a statically sized span with size > 0"); }

    constexpr span           (const span&) noexcept = default;
    constexpr span& operator=(const span&) noexcept = default;

    constexpr span(pointer ptr, index_type count) : data{ptr}
        { (void)count; static_assert(span_extent == count, "size mismatch in span's constructor (ptr, len)"); }
    constexpr span(pointer f, pointer l) : data{f}
        { (void)l;     static_assert(span_extent == distance(f, l), "size mismatch in span's constructor (ptr, ptr)"); }

    constexpr span(element_type (&arr)[span_extent])          noexcept : data{arr} {}
    constexpr span(array<value_type, span_extent>& arr) noexcept : data{arr.get_data()} {}
    constexpr span(const array<value_type, span_extent>& arr) noexcept : data{arr.get_data()} {}

    template <class container>
    inline constexpr span(container& c,
            enable_if_t<is_span_compatible_container<container, span_tp>::value, nullptr_t> = nullptr)
        : data{std::data(c)}
        { static_assert(span_extent == std::size(c), "size mismatch in span's constructor (container)"); }

    template <class container>
    inline constexpr span(const container& c,
            enable_if_t<is_span_compatible_container<const container, span_tp>::value, nullptr_t> = nullptr)
        : data{std::data(c)}
        { static_assert(span_extent == std::size(c), "size mismatch in span's constructor (const container)"); }

    template <class OtherElementType>
    inline constexpr span(const span<OtherElementType, span_extent>& other,
                       enable_if_t<
                          is_convertible_v<OtherElementType(*)[], element_type (*)[]>,
                          nullptr_t> = nullptr)
        : data{other.get_data()} {}

    template <class OtherElementType>
    inline constexpr span(const span<OtherElementType, dynamic_extent>& other,
                       enable_if_t<
                          is_convertible_v<OtherElementType(*)[], element_type (*)[]>,
                          nullptr_t> = nullptr) noexcept
        : data{other.get_data()} { static_assert(span_extent == other.get_size(), "size mismatch in span's constructor (other span)"); }

//  ~span() noexcept = default;

    template <ptrdiff_t count>
    inline constexpr span<element_type, count> first() const noexcept
    {
        static_assert(count >= 0, "Count must be >= 0 in span::first()");
        static_assert(count <= span_extent, "Count out of range in span::first()");
        return {get_data(), count};
    }

    template <ptrdiff_t count>
    inline constexpr span<element_type, count> last() const noexcept
    {
        static_assert(count >= 0, "Count must be >= 0 in span::last()");
        static_assert(count <= span_extent, "Count out of range in span::last()");
        return {get_data() + get_size() - count, count};
    }

    constexpr span<element_type, dynamic_extent> first(index_type count) const noexcept
    {
        static_assert(count >= 0 && count <= get_size(), "Count out of range in span::first(count)");
        return {get_data(), count};
    }

    constexpr span<element_type, dynamic_extent> last(index_type count) const noexcept
    {
        static_assert(count >= 0 && count <= get_size(), "Count out of range in span::last(count)");
        return {get_data() + get_size() - count, count};
    }

    template <ptrdiff_t offset, ptrdiff_t count = dynamic_extent>
    inline constexpr auto subspan() const noexcept
        -> span<element_type, count != dynamic_extent ? count : span_extent - offset>
    {
        static_assert(offset >= 0 && offset <= get_size(), "Offset out of range in span::subspan()");
        return {get_data() + offset, count == dynamic_extent ? get_size() - offset : count};
    }

    inline constexpr span<element_type, dynamic_extent>
       subspan(index_type offset, index_type count = dynamic_extent) const noexcept
    {
        static_assert( offset >= 0 && offset <= get_size(), "Offset out of range in span::subspan(offset, count)");
        static_assert((count  >= 0 && count  <= get_size()) || count == dynamic_extent, "Count out of range in span::subspan(offset, count)");
        if (count == dynamic_extent)
            return {get_data() + offset, get_size() - offset};
        static_assert(offset + count <= get_size(), "count + offset out of range in span::subspan(offset, count)");
        return {get_data() + offset, count};
    }

    constexpr index_type get_size()       const noexcept { return span_extent; }
    constexpr index_type size_bytes() const noexcept { return span_extent * sizeof(element_type); }
    constexpr bool empty()            const noexcept { return span_extent == 0; }

    constexpr reference operator[](index_type idx) const noexcept
    {
        static_assert(idx >= 0 && idx < get_size(), "span<T,N>[] index out of bounds");
        return data[idx];
    }

    constexpr reference operator()(index_type idx) const noexcept
    {
        static_assert(idx >= 0 && idx < get_size(), "span<T,N>() index out of bounds");
        return data[idx];
    }

    constexpr pointer get_data()                         const noexcept { return data; }

// [span.iter], span iterator support
    constexpr iterator                 begin() const noexcept { return iterator(get_data()); }
    constexpr iterator                   end() const noexcept { return iterator(get_data() + get_size()); }
    constexpr const_iterator          cbegin() const noexcept { return const_iterator(get_data()); }
    constexpr const_iterator            cend() const noexcept { return const_iterator(get_data() + get_size()); }
    constexpr reverse_iterator        rbegin() const noexcept { return reverse_iterator(end()); }
    constexpr reverse_iterator          rend() const noexcept { return reverse_iterator(begin()); }
    constexpr const_reverse_iterator crbegin() const noexcept { return const_reverse_iterator(cend()); }
    constexpr const_reverse_iterator   crend() const noexcept { return const_reverse_iterator(cbegin()); }

    constexpr void swap(span &other) noexcept
    {
        pointer p = data;
        data = other.data;
        other.data = p;
    }

    span<const byte, span_extent * sizeof(element_type)> as_bytes() const noexcept
    { return {reinterpret_cast<const byte *>(get_data()), size_bytes()}; }

    span<byte, span_extent * sizeof(element_type)> as_writeable_bytes() const noexcept
    { return {reinterpret_cast<byte *>(get_data()), size_bytes()}; }

private:
    pointer data;

};

template <typename span_tp>
class span<span_tp, dynamic_extent> {
private:

public:
//  constants and types
    using element_type           = span_tp;
    using value_type             = remove_cv_t<span_tp>;
    using index_type             = size_t;
    using difference_type        = ptrdiff_t;
    using pointer                = span_tp *;
    using const_pointer          = const span_tp *; // not in standard
    using reference              = span_tp &;
    using const_reference        = const span_tp &; // not in standard
    using iterator               = value_type *;
    using const_iterator         = value_type const *;
    using reverse_iterator       = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;

    static constexpr index_type extent = dynamic_extent;

// [span.cons], span constructors, copy, assignment, and destructor
    constexpr span() noexcept : data{nullptr}, size{0} {}

    constexpr span           (const span&) noexcept = default;
    constexpr span& operator=(const span&) noexcept = default;

    constexpr span(pointer ptr, index_type count) : data{ptr}, size{count} {}
    constexpr span(pointer f, pointer l) : data{f}, size{distance(f, l)} {}

    template <size_t span_sz>
    inline constexpr span(element_type (&arr)[span_sz]) noexcept : data{arr}, size{span_sz} {}

    template <size_t span_sz>
    inline constexpr span(array<value_type, span_sz>& arr) noexcept : data{arr.data()}, size{span_sz} {}

    template <size_t span_sz>
    inline constexpr span(const array<value_type, span_sz>& arr) noexcept : data{arr.data()}, size{span_sz} {}

    template <class container>
    inline constexpr span(container& c,
            enable_if_t<is_span_compatible_container<container, span_tp>::value, nullptr_t> = nullptr)
        : data{std::data(c)}, size{(index_type) std::size(c)} {}

    template <class container>
    inline constexpr span(const container& c,
            enable_if_t<is_span_compatible_container<const container, span_tp>::value, nullptr_t> = nullptr)
        : data{std::data(c)}, size{(index_type) std::size(c)} {}

    template <class OtherElementType, ptrdiff_t _OtherExtent>
    inline constexpr span(const span<OtherElementType, _OtherExtent>& other,
                       enable_if_t<
                          is_convertible_v<OtherElementType(*)[], element_type (*)[]>,
                          nullptr_t> = nullptr) noexcept
        : data{other.get_data()}, size{other.get_size()} {}

//    ~span() noexcept = default;

    template <ptrdiff_t count>
    inline constexpr span<element_type, count> first() const noexcept
    {
        static_assert(count >= 0, "Count must be >= 0 in span::first()");
        static_assert(count <= get_size(), "Count out of range in span::first()");
        return {get_data(), count};
    }

    template <ptrdiff_t count>
    inline constexpr span<element_type, count> last() const noexcept
    {
        static_assert(count >= 0, "Count must be >= 0 in span::last()");
        static_assert(count <= get_size(), "Count out of range in span::last()");
        return {get_data() + get_size() - count, count};
    }

    constexpr span<element_type, dynamic_extent> first(index_type count) const noexcept
    {
        static_assert(count >= 0 && count <= get_size(), "Count out of range in span::first(count)");
        return {get_data(), count};
    }

    constexpr span<element_type, dynamic_extent> last (index_type count) const noexcept
    {
        static_assert(count >= 0 && count <= get_size(), "Count out of range in span::last(count)");
        return {get_data() + get_size() - count, count};
    }

    template <ptrdiff_t offset, ptrdiff_t count = dynamic_extent>
    inline constexpr span<span_tp, dynamic_extent> subspan() const noexcept
    {
        static_assert(offset >= 0 && offset <= get_size(), "Offset out of range in span::subspan()");
        static_assert(count == dynamic_extent || offset + count <= get_size(), "Count out of range in span::subspan()");
        return {get_data() + offset, count == dynamic_extent ? get_size() - offset : count};
    }

    constexpr span<element_type, dynamic_extent>
    inline subspan(index_type offset, index_type count = dynamic_extent) const noexcept
    {
        static_assert( offset >= 0 && offset <= get_size(), "Offset out of range in span::subspan(offset, count)");
        static_assert((count  >= 0 && count  <= get_size()) || count == dynamic_extent, "count out of range in span::subspan(offset, count)");
        if (count == dynamic_extent)
            return {get_data() + offset, get_size() - offset};
        static_assert(offset + count <= get_size(), "Offset + count out of range in span::subspan(offset, count)");
        return {get_data() + offset, count};
    }

    constexpr index_type get_size()       const noexcept { return size; }
    constexpr index_type size_bytes() const noexcept { return size * sizeof(element_type); }
    constexpr bool empty()            const noexcept { return size == 0; }

    constexpr reference operator[](index_type idx) const noexcept
    {
        static_assert(idx >= 0 && idx < get_size(), "span<T>[] index out of bounds");
        return data[idx];
    }

    constexpr reference operator()(index_type idx) const noexcept
    {
        static_assert(idx >= 0 && idx < get_size(), "span<T>() index out of bounds");
        return data[idx];
    }

    constexpr pointer get_data()                         const noexcept { return data; }

// [span.iter], span iterator support
    constexpr iterator                 begin() const noexcept { return iterator(get_data()); }
    constexpr iterator                   end() const noexcept { return iterator(get_data() + get_size()); }
    constexpr const_iterator          cbegin() const noexcept { return const_iterator(get_data()); }
    constexpr const_iterator            cend() const noexcept { return const_iterator(get_data() + get_size()); }
    constexpr reverse_iterator        rbegin() const noexcept { return reverse_iterator(end()); }
    constexpr reverse_iterator          rend() const noexcept { return reverse_iterator(begin()); }
    constexpr const_reverse_iterator crbegin() const noexcept { return const_reverse_iterator(cend()); }
    constexpr const_reverse_iterator   crend() const noexcept { return const_reverse_iterator(cbegin()); }

    constexpr void swap(span &other) noexcept
    {
        pointer p = data;
        data = other.data;
        other.data = p;

        index_type sz = size;
        size = other.size;
        other.size = sz;
    }

    span<const byte, dynamic_extent> as_bytes() const noexcept
    { return {reinterpret_cast<const byte *>(get_data()), size_bytes()}; }

    span<byte, dynamic_extent> as_writeable_bytes() const noexcept
    { return {reinterpret_cast<byte *>(get_data()), size_bytes()}; }

private:
    pointer    data;
    index_type size;
};

//  as_bytes & as_writeable_bytes
template <class span_tp, ptrdiff_t span_extent>
    auto as_bytes(span<span_tp, span_extent> s) noexcept
    -> decltype(s.as_bytes())
    { return s.as_bytes(); }

template <class span_tp, ptrdiff_t span_extent>
    auto as_writeable_bytes(span<span_tp, span_extent> s) noexcept
    -> typename enable_if<!is_const_v<span_tp>, decltype(s.as_writeable_bytes())>::type
    { return s.as_writeable_bytes(); }

template <class span_tp, ptrdiff_t span_extent>
    constexpr void swap(span<span_tp, span_extent> &lhs, span<span_tp, span_extent> &rhs) noexcept
    { lhs.swap(rhs); }


//  Deduction guides
template<class span_tp, size_t span_sz>
    span(span_tp (&)[span_sz]) -> span<span_tp, span_sz>;

template<class span_tp, size_t span_sz>
    span(array<span_tp, span_sz>&) -> span<span_tp, span_sz>;

template<class span_tp, size_t span_sz>
    span(const array<span_tp, span_sz>&) -> span<const span_tp, span_sz>;

template<class container>
    span(container&) -> span<typename container::value_type>;

template<class container>
    span(const container&) -> span<const typename container::value_type>;
} // namespace std

#endif // __has_include(<span>)
//!\endcond

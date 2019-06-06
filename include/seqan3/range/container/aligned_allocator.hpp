// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 * \brief Provides seqan3::aligned_allocator.
 */

#pragma once

#include <limits>
#include <memory>
#include <type_traits>

#include <seqan3/core/platform.hpp>

namespace seqan3
{

/*!\brief Allocates uninitialized storage whose memory-alignment is specified by *alignment*.
 * \tparam value_t     \copydoc aligned_allocator::value_type
 * \tparam alignment_v \copydoc aligned_allocator::alignment
 * \ingroup container
 *
 * \details
 *
 * \include test/snippet/range/container/aligned_allocator.cpp
 *
 * Will output something like:
 * ```console
 * Item: 1 (0x55d5d0722f00, 128-byte aligned offset: 0)
 * Item: 2 (0x55d5d0722f02, 128-byte aligned offset: 2)
 * Item: 3 (0x55d5d0722f04, 128-byte aligned offset: 4)
 * Item: 4 (0x55d5d0722f06, 128-byte aligned offset: 6)
 * Item: 5 (0x55d5d0722f08, 128-byte aligned offset: 8)
 * Item: 1 (0x55d5d0722f40, unaligned start: 64)
 * Item: 2 (0x55d5d0722f42, unaligned start: 66)
 * Item: 3 (0x55d5d0722f44, unaligned start: 68)
 * Item: 4 (0x55d5d0722f46, unaligned start: 70)
 * Item: 5 (0x55d5d0722f48, unaligned start: 72)
 * Item: 1 (0x55d5d0723000, 256-byte aligned offset: 0)
 * Item: 2 (0x55d5d0723004, 256-byte aligned offset: 4)
 * Item: 3 (0x55d5d0723008, 256-byte aligned offset: 8)
 * Item: 4 (0x55d5d072300c, 256-byte aligned offset: 12)
 * Item: 5 (0x55d5d0723010, 256-byte aligned offset: 16)
 * ```
 *
 * As you can see, in the case of the aligned_allocator it is guaranteed that the
 * first element in the vector starts at offset 0.
 *
 * \see http://en.cppreference.com/w/cpp/concept/Allocator
 * \see http://en.cppreference.com/w/cpp/memory/c/aligned_alloc
 */
template <typename value_t, size_t alignment_v = __STDCPP_DEFAULT_NEW_ALIGNMENT__>
class aligned_allocator
{
public:
    //!\brief The memory-alignment of the allocation
    static constexpr size_t alignment = alignment_v;

    //!\brief The value type of the allocation
    using value_type = value_t;
    //!\brief The pointer type of the allocation
    using pointer = value_type*;
    //!\brief The difference type of the allocation
    using difference_type = typename std::pointer_traits<pointer>::difference_type;
    //!\brief The size type of the allocation
    using size_type = std::make_unsigned_t<difference_type>;

    //!\brief Are any two allocators of the same aligned_allocator type always compare equal?
    using is_always_equal = std::true_type;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    aligned_allocator()                                     = default; //!< Defaulted.
    aligned_allocator(aligned_allocator const &)            = default; //!< Defaulted.
    aligned_allocator(aligned_allocator &&)                 = default; //!< Defaulted.
    aligned_allocator& operator=(aligned_allocator const &) = default; //!< Defaulted.
    aligned_allocator& operator=(aligned_allocator &&)      = default; //!< Defaulted.
    ~aligned_allocator()                                    = default; //!< Defaulted.

    //!\brief Copy constructor with different value type.
    template <class other_value_type>
    constexpr aligned_allocator(aligned_allocator<other_value_type, alignment> const &) noexcept
    {}
    //!\}

    /*!\brief Allocates `n * sizeof(T)` bytes of uninitialized storage by calling std::aligned_alloc, but it is
     * unspecified when and how this function is called.
     * \throws Throws std::bad_alloc if allocation fails.
     * \sa https://en.cppreference.com/w/cpp/memory/allocator/allocate
     */
    [[nodiscard]]
    pointer allocate(size_type n)
    {
        constexpr size_type max_size = std::numeric_limits<size_type>::max() / sizeof(value_type);
        if (n > max_size)
            throw std::bad_alloc();

        // NOTE: On macOS glibc does not implement aligned_alloc, so we need to fallback to posix_memalign instead.
#if defined(__APPLE__) && (!defined(_GLIBCXX_HAVE_ALIGNED_ALLOC) && !defined(_ISOC11_SOURCE))
        void * p{};
        if (int res = posix_memalign(&p, alignment, n * sizeof(value_type)); res == 0 && p != nullptr)
            return static_cast<pointer>(p);
#else
        // NOTE:
        // Allocate size bytes of uninitialized storage whose alignment is
        // specified by alignment. The size parameter must be an integral
        // multiple of alignment.
        //
        // Passing a size which is not an integral multiple of alignment or an
        // alignment which is not valid or not supported by the implementation
        // causes the function to fail and return a null pointer (C11, as
        // published, specified undefined behavior in this case, this was
        // corrected by DR 460).
        // http://en.cppreference.com/w/cpp/memory/c/aligned_alloc
        if (auto p = static_cast<pointer>(std::aligned_alloc(alignment, n * sizeof(value_type))))
            return p;
#endif

        throw std::bad_alloc();
    }

    /*!\brief Deallocates the storage referenced by the pointer p, which must be a pointer obtained by an earlier call
     * to allocate().
     * \details
     *
     * The argument n must be equal to the first argument of the call to allocate() that originally produced p;
     * otherwise, the behavior is undefined.
     *
     * Calls std::free, but it is unspecified when and how it is called.
     * \sa https://en.cppreference.com/w/cpp/memory/allocator/deallocate
     */
    void deallocate(pointer p, size_type) noexcept
    {
        std::free(p);
    }

    /*!\brief The aligned_allocator member template class aligned_allocator::rebind provides a way to obtain an
     * allocator for a different type.
     * \tparam new_value_type The other value type.
     */
    template<typename new_value_type>
    struct rebind
    {
        //!\brief The type of the allocator for a different value type.
        using other = aligned_allocator<new_value_type, alignment>;
    };

    /*!\name Comparison operators
     * \{
     */
    //!\brief Returns true if the memory-alignment matches.
    template <class value_type2, size_t alignment2>
    constexpr bool operator==(aligned_allocator<value_type2, alignment2> const &) noexcept
    {
        return alignment == alignment2;
    }

    //!\brief Returns false if the memory-alignment mismatches.
    template <class value_type2, size_t alignment2>
    constexpr bool operator!=(aligned_allocator<value_type2, alignment2> const &) noexcept
    {
        return alignment != alignment2;
    }
    //!\}
};

} // namespace seqan3

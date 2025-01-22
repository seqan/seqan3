// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

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
 * \tparam value_t     The value type of the allocation.
 * \tparam alignment_v The memory-alignment of the allocation; defaults to `__STDCPP_DEFAULT_NEW_ALIGNMENT__`.
 * \ingroup utility_container
 *
 * \details
 *
 * This class allocates memory at the given `alignment_v` offset. This makes sure that the allocated memory
 * starts at a memory offset equal to some multiple of the word size. More formally, a memory address `a`, is said to
 * be `n`-byte aligned when `n` is a power of two and `a` is a multiple of `n` bytes.
 *
 * If the specified `alignment` is not supported (e.g. alignments that are not a power of two) by the
 * used allocation method a std::bad_alloc exception will be thrown. For requested alignments larger than
 * `__STDCPP_DEFAULT_NEW_ALIGNMENT__`, also called new-extended alignments, the storage will have the alignment
 * specified by the value `alignment`. Otherwise, the storage is aligned for any object that does not have new-extended
 * alignment, e.g. `int` or `double`, and is of the requested size.
 *
 * \include test/snippet/utility/container/aligned_allocator.cpp
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
 * \see https://en.cppreference.com/w/cpp/named_req/Allocator
 * \see https://en.cppreference.com/w/cpp/memory/c/aligned_alloc
 *
 * \noapi{Utility for allocator-aware container.}
 */
template <typename value_t, size_t alignment_v = __STDCPP_DEFAULT_NEW_ALIGNMENT__>
class aligned_allocator
{
public:
    //!\brief The memory-alignment of the allocation.
    static constexpr size_t alignment = alignment_v;

    //!\brief The value type of the allocation.
    using value_type = value_t;
    //!\brief The pointer type of the allocation.
    using pointer = value_type *;
    //!\brief The difference type of the allocation.
    using difference_type = typename std::pointer_traits<pointer>::difference_type;
    //!\brief The size type of the allocation.
    using size_type = std::make_unsigned_t<difference_type>;

    //!\brief Do any two allocators of the same aligned_allocator type always compare equal?
    using is_always_equal = std::true_type;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    aligned_allocator() = default;                                      //!< Defaulted.
    aligned_allocator(aligned_allocator const &) = default;             //!< Defaulted.
    aligned_allocator(aligned_allocator &&) = default;                  //!< Defaulted.
    aligned_allocator & operator=(aligned_allocator const &) = default; //!< Defaulted.
    aligned_allocator & operator=(aligned_allocator &&) = default;      //!< Defaulted.
    ~aligned_allocator() = default;                                     //!< Defaulted.

    //!\brief Copy constructor with different value type and alignment.
    template <class other_value_type, size_t other_alignment>
    constexpr aligned_allocator(aligned_allocator<other_value_type, other_alignment> const &) noexcept
    {}
    //!\}

    /*!\brief Allocates sufficiently large memory to hold `n` many elements of `value_type`.
     *
     * \param[in] n The number of elements for which to allocate the memory.
     *
     * \returns The pointer to the first block of allocated memory.
     *
     * \throws Throws std::bad_alloc if allocation fails, i.e. either the call to the throwing version of
     *         operator new function throws or the requested memory exceeds the maximal number of elements to allocate.
     *
     * \details
     *
     * Allocates `n * sizeof(value_type)` bytes of uninitialized storage by calling
     * [operator new](https://en.cppreference.com/w/cpp/memory/new/operator_new).
     * If the given `alignment` is bigger than
     * [__STDCPP_DEFAULT_NEW_ALIGNMENT__](https://en.cppreference.com/w/cpp/memory/new/align_val_t), the alignment
     * aware operator new that takes as second argument the desired alignment of type std::align_val_t is used.
     *
     * \note We call the new operator with the semantic requirements that the c++ standard specifies/demands, but be
     *       aware that users can overload any (global) `operator new` that might not adhere to the standard and might
     *       cause std::bad_alloc or unaligned pointers.
     *
     * \sa https://en.cppreference.com/w/cpp/memory/allocator/allocate
     *
     * ### Thread safety
     *
     * Thread-safe.
     *
     * ### Exception
     *
     * Strong exception guarantee.
     */
    [[nodiscard]] pointer allocate(size_type const n) const
    {
        constexpr size_type max_size = std::numeric_limits<size_type>::max() / sizeof(value_type);
        if (n > max_size)
            throw std::bad_alloc{};

        size_t bytes_to_allocate = n * sizeof(value_type);
        if constexpr (alignment <= __STDCPP_DEFAULT_NEW_ALIGNMENT__)
            return static_cast<pointer>(::operator new(bytes_to_allocate));
        else // Use alignment aware allocator function.
            return static_cast<pointer>(::operator new(bytes_to_allocate, static_cast<std::align_val_t>(alignment)));
    }

    /*!\brief Deallocates the storage referenced by the pointer p, which must be a pointer obtained by an earlier call
     * to seqan3::aligned_allocator::allocate.
     *
     * \param[in] p The pointer to the memory to be deallocated.
     * \param[in] n The number of elements to be deallocated.
     *
     * \details
     *
     * The argument `n` must be equal to the first argument of the call to seqan3::aligned_allocator::allocate that
     * originally produced `p`, otherwise the behavior is undefined. This function calls
     * [operator delete](https://en.cppreference.com/w/cpp/memory/new/operator_delete) to deallocate the memory of
     * specified size.
     * If the given `alignment` is bigger than
     * [__STDCPP_DEFAULT_NEW_ALIGNMENT__](https://en.cppreference.com/w/cpp/memory/new/align_val_t) the alignment
     * aware operator delete that takes as third argument the alignment as std::align_val_t.
     *
     * \sa https://en.cppreference.com/w/cpp/memory/allocator/deallocate
     *
     * ### Thread safety
     *
     * Thread-safe.
     *
     * ### Exception
     *
     * Nothrow guarantee.
     */
    void deallocate(pointer const p, size_type const n) const noexcept
    {
        size_t bytes_to_deallocate = n * sizeof(value_type);

        if constexpr (alignment <= __STDCPP_DEFAULT_NEW_ALIGNMENT__)
            ::operator delete(p, bytes_to_deallocate);
        else // Use alignment aware deallocator function.
            ::operator delete(p, bytes_to_deallocate, static_cast<std::align_val_t>(alignment));
    }

    /*!\brief The aligned_allocator member template class aligned_allocator::rebind provides a way to obtain an
     *        allocator for a different type.
     *
     * \tparam new_value_type The other value type.
     *
     * \details
     *
     * If the alignment of the new type exceeds the alignment of the current allocator, the larger alignment will
     * be used.
     */
    template <typename new_value_type>
    struct rebind
    {
        //!\brief The alignment for the rebound allocator.
        static constexpr size_t other_alignment = std::max(alignof(new_value_type), alignment);
        //!\brief The type of the allocator for a different value type.
        using other = aligned_allocator<new_value_type, other_alignment>;
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

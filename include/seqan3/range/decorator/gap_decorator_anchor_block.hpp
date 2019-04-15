// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::gap_decorator_anchor_block.
 * \author Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
 */

#pragma once

#include <numeric>
#include <tuple>
#include <type_traits>
#include <vector>

#include <seqan3/alignment/exception.hpp>
#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/gap/gap.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/range/container/concept.hpp>
#include <seqan3/range/detail/random_access_iterator.hpp>
#include <seqan3/std/ranges>

namespace seqan3
{

/*!\brief A gap decorator allows the annotation of sequences with gap symbols
 *        while leaving the underlying sequence unmodified.
 * \tparam inner_type The type of range that will be decorated with gaps; must model std::ranges::RandomAccessRange
 *                    and std::ranges::SizedRange.
 * \implements seqan3::aligned_sequence_concept
 * \ingroup decorator
 *
 * \details
 *
 * This class may be used whenever you want to store or compute an alignment. The
 * underlying (ungapped) sequence remains unmodified, and is augmented with gap
 * information. The seqan3::gap_decorator_anchor_block behaves just like a random
 * access container of over a gapped alphabet when iterating over it,
 * inserting/erasing gaps or accessing a position. The only difference lies in the
 * performance and size overhead (see below).
 *
 * ### Performance
 *
 * **n** The length of the underlying sequence.
 * **k** The number of contiguous gaps (not gap symbols).
 * **b** The block width.
 *
 * |            | access next | random access               | gap insert/erase at end      | gap insert/erase random           | size overhead            |
 * |------------|-------------|---------------------------- |------------------------------|-----------------------------------|--------------------------|
 * | decorator  | \f$O(1)\f$  | \f$O(\frac{n}{b}\log(b))\f$ | \f$O(\frac{n}{b}\log(b))\f$  | \f$O(\frac{n}{b}(\log(b) + 1))\f$ | \f$O(\frac{n}{b} + k)\f$ |
 * | vector     | \f$O(1)\f$  | \f$O(1)\f$                  | \f$O(1)\f$                   | \f$O(n)\f$                        | \f$O(n)\f$               |
 *
 * The *size overhead* refers to the space that is needed when using each of the data structures in addition to an
 * already existing ungapped sequence.
 *
 * ### Implementation details
 *
 * This decorator represents gaps as anchors by position and gap length. The
 * position is relative to the underlying, ungapped sequence which is given by
 * reference. The list of anchor gaps is sorted and hierarchically organized into
 * blocks of fixed maximum length (set via template parameter `block_size`). To
 * aid faster search for a given virtual address (`pos`), accumulated gap lengths
 * are stored in an extra vector in the size of the number of blocks (plus one
 * for tailing gaps). Random access (needed for reading, gap insertion/erasure)
 * is always performed in two phases: firstly, binary search on the accumulative
 * statistics (`gap_sums`) to identify the block holding the anchor gap following
 * or surrounding the virtual address, secondly, linear search within the block
 * to identify the lower bounding anchor gap.
 * Upon insertion of gap symbols either a new gap is inserted into its designated
 * block or an existing one extended. Finally, the tail of the vector of accumulated
 * gap sums has to be updated including the index of the modified block.
 * Gap erasure works analogously.
 *
 * ### The seqan3::gap_decorator_anchor_block::iterator type
 *
 * \attention The iterator of the seqan3::gap_decorator_anchor_block does not model the
 *            [Cpp17InputIterator](https://en.cppreference.com/w/cpp/named_req/InputIterator) requirements of the
 *            STL because dereferencing the iterator returns a proxy and no operator-> is provided.
 *            Note that it does model the std::ranges::InputIterator.
 *
 */
template <std::ranges::ViewableRange inner_type, unsigned int block_size=32>
//!\cond
    requires std::ranges::RandomAccessRange<inner_type> && std::ranges::SizedRange<inner_type> &&
             (std::is_const_v<std::remove_reference_t<inner_type>> || std::ranges::View<inner_type>)
//!\endcond
class gap_decorator_anchor_block
{
private:
    /*!\brief The iterator that moves over the seqan3::gap_decorator_anchor_block.
     *
     * \details
     *
     * This iterator returns values when dereferenced, not references, i.e. it does
     * not satisfy the semantic requirements of [LegacyForwardIterator](https://en.cppreference.com/w/cpp/named_req/ForwardIterator).
     * It does model the C++20 std::BidirectionalIterator (and std::ForwardIterator
     * implicitly).
     */
    class gap_decorator_anchor_block_iterator
    {
    private:
        friend class gap_decorator_anchor_block;

        //!\brief Pointer to the underlying container structure.
        typename std::add_pointer_t<gap_decorator_anchor_block const> host{nullptr};
        //!\brief Stores the virtual position index for the seqan3::gap_decorator_anchor_block.
        typename gap_decorator_anchor_block::size_type pos{0u};
        //!\brief Stores the physical position in the ungapped/underlying view.
        int64_t ungapped_view_pos{0}; // must be signed because we need this value to be -1 in case of leading gaps.
        //!\brief Stores current gap block index designated for `pos` in the alignment space.
        typename gap_decorator_anchor_block::size_type block_id{0};
        //!\brief Stores the lower bounding gap index within the block for `pos`.
        //!       A lower bounding gap is either enclosing, succeeding (if there
        //!       is no enclosing one) or points to the end of the block (if there
        //!       is no succeeding gap).
        typename gap_decorator_anchor_block::size_type gap_id{0};
        //!\brief Caches whether the iterator points to a gap (true) or not (false).
        //!     In case the current position is a gap symbol, `gap_id` points to the
        //!     corresponding gap, else `gap_id` points to the succeeding gap in the
        //!     designated block or to its end if there is no gap.
        bool is_at_gap{true};
        //!\brief An accumulator for all preceeding (and enclosing gaps in case is_at_gap is true).
        typename gap_decorator_anchor_block::size_type gap_acc{0};

        //!\brief A helper function that performs the random access into the gap
        //! block list and updates all member variables.
        // TODO: set ungapped_view_pos
        void jump(typename gap_decorator_anchor_block::size_type const new_pos)
        {
            assert(new_pos <= host->size());
            pos = new_pos;
            // computes the right bounding virtual position of block given by its index
            std::function<size_type (size_type)> right_bound =
                [=](size_type i) { return (i+1)*block_size +
                    host->gap_sums[i]; };

            size_type left = 0, right = host->gap_sums.size(), mid = host->gap_sums.size()/2;
            // number of gaps and alphabet symbols until end of this block
            size_type num_symbols = right_bound(mid);
            // number of gaps and alphabet symbols until end of preceeding block
            size_type num_symbols_pred = (mid) ? right_bound(mid-1) : 0;

            // binary search to identify designated block
            while (pos < num_symbols_pred || pos >= num_symbols)
            {
                if (pos < num_symbols_pred && mid)
                    right = mid; //mid /= 2;
                else //if (mid < gap_sums.size() - 1)
                    left = mid; //mid = mid + (gap_sums.size() - mid)/2;  // TODO: check for floor for fractional numbers
                mid = (left + right)/2;
                num_symbols = right_bound(mid);
                num_symbols_pred = (mid) ? right_bound(mid-1) : 0;
            }
            // locate upper bounding gap within block, if there is no upper bound gap_id points to end
            gap_id = 0;
            gap_acc = (mid) ? host->gap_sums[mid-1] : 0;
            // accumulate gaps as long 'pos' refers to position before or inside the current gap
            while (gap_id < host->gap_block_list[mid].size() && pos >= gap_acc + host->gap_block_list[mid][gap_id].first)
            {
                gap_acc += host->gap_block_list[mid][gap_id].second;
                ++gap_id;
            }
            block_id = mid;
            is_at_gap = false;
            if (!gap_id && gap_id < host->gap_block_list[block_id].size() &&
                    pos >= gap_acc + host->gap_block_list[block_id][gap_id].first)
                is_at_gap = true;
            // increment gap_id if pos refers to position inside gap, i.e.
            if (gap_id && pos >= gap_acc + host->gap_block_list[mid][gap_id-1].first -
                    host->gap_block_list[mid][gap_id-1].second && pos < gap_acc + host->gap_block_list[mid][gap_id-1].first)
            {
                --gap_id;
                is_at_gap = true;
            }
            ungapped_view_pos = pos - gap_acc;
        }

    public:
        /*!\name Member types
         * \brief Make the parent's member types visible.
         * \{
         */
        //!\brief Type for distances between iterators.
        using difference_type = typename gap_decorator_anchor_block::difference_type;
        //!\brief Value type of container elements.
        using value_type = typename gap_decorator_anchor_block::value_type;
        //!\brief Const-reference type defined by container (which equals the reference type).
        using reference = typename gap_decorator_anchor_block::const_reference; // = reference
        //!\brief Equals reference type.
        using const_reference = reference;
        //!\brief Pointer of container value type.
        using pointer = value_type *;
        //!\brief Tag this class as a bidirectional iterator since random access is possible but not in constant time.
        using iterator_category = std::bidirectional_iterator_tag;
        //!\}

        /*!\name Constructors/Destructors
         * \{
         */
        //!\brief Default constructor.
        constexpr gap_decorator_anchor_block_iterator() = default;
        //!\brief Copy constructor.
        constexpr gap_decorator_anchor_block_iterator(gap_decorator_anchor_block_iterator const &) = default;
        //!\brief Copy construction via assignment.
        constexpr gap_decorator_anchor_block_iterator & operator=(gap_decorator_anchor_block_iterator const &) = default;
        //!\brief Move constructor.
        constexpr gap_decorator_anchor_block_iterator (gap_decorator_anchor_block_iterator &&) = default;
        //!\brief Move assignment.
        constexpr gap_decorator_anchor_block_iterator & operator=(gap_decorator_anchor_block_iterator &&) = default;
        //!\brief Use default deconstructor.
        ~gap_decorator_anchor_block_iterator() = default;

        //!\brief Construct from seqan3::gap_decorator_anchor_block and initialise members.
        explicit constexpr gap_decorator_anchor_block_iterator(gap_decorator_anchor_block const & host_) noexcept :
            host(&host_), pos{0}, ungapped_view_pos{0}, block_id{0}, gap_id{0}, gap_acc{0}
        {
            is_at_gap = false;
            if (host_.gap_block_list.size() && host_.gap_block_list[0].size() && !(*host_.gap_block_list[0].begin()).first)
            {
                is_at_gap = true;
                gap_acc = (*host_.gap_block_list[0].begin()).second;
                ungapped_view_pos = -1;
            }
        }

        //!\brief Construct from seqan3::gap_decorator_anchor_block and explicit position.
        constexpr gap_decorator_anchor_block_iterator(gap_decorator_anchor_block const & host_,
                                                    typename gap_decorator_anchor_block::size_type const pos_) noexcept :
             host(&host_)
        {
            jump(pos_); // random access to pos
        }
        //!\}

        /*!\name Arithmetic operators
         * \{
        */
        //!\brief Pre-increment, returns updated iterator.
        constexpr gap_decorator_anchor_block_iterator & operator++() noexcept
        {
            assert(host); // host is set
            ++pos;
            if (is_at_gap)
            {
                // rightmost gap position of enclosing anchor gap
                size_type gap_end = gap_acc + host->gap_block_list[block_id][gap_id].first;
                if (pos < gap_end)
                {
                    // is_at_gap and indices remain unmodified
                }
                else
                {
                    is_at_gap = false;
                    // exploit invariant that gap ending in a block is always followed virtually by non-gap symbol
                    ++gap_id;  // points now to succeeding gap or end of block
                }
            }
            else // !is_at_gap
            {
                // ungapped view to block mapping
                size_type p_pos = pos - gap_acc - block_size*block_id;
                // switch to next block
                if (!(p_pos % block_size))
                {
                    ++block_id;
                    gap_id = 0;
                    if (host->gap_block_list[block_id].size() && !((*host->gap_block_list[block_id].begin()).first % block_size))
                        is_at_gap = true;
                }
                else // stay in same block
                {
                    if (!host->gap_block_list[block_id].size()
                        || gap_id == host->gap_block_list[block_id].size())
                    {
                        // no succeeding gap, do nothing
                    }
                    else
                    {
                        // succeeding gap start position
                        size_type gap_start = gap_acc + host->gap_block_list[block_id][gap_id].first;
                        // pos points now into gap start
                        if (gap_start <= pos) //==,  succeeding gap coincides with updated position
                        {
                            is_at_gap = true;
                            gap_acc += host->gap_block_list[block_id][gap_id].second;
                        }
                    }
                }
            }
            ungapped_view_pos = pos - gap_acc;
            return *this;
        }

        //!\brief Pre-decrement, returns updated iterator.
        constexpr gap_decorator_anchor_block_iterator & operator--() noexcept
        {
            assert(host); // host is set
            assert(pos); // pos=0 not decrementable
            --pos;
            // is at gap, but may now point to ungapped view letter preceeding current gap
            if (is_at_gap)
            {
                size_type current_gap_start = host->gap_block_list[block_id][gap_id].first + gap_acc - host->gap_block_list[block_id][gap_id].second;
                // drop gap
                if (current_gap_start > pos)
                {
                    // substract formerly enclosing gap length
                    gap_acc -= host->gap_block_list[block_id][gap_id].second;
                    is_at_gap = false;
                    // switch pointers to now end of preceeding block
                    if (!gap_id)
                    {
                        --block_id; // note: block_id > 0, otherwise original pos = 0
                        gap_id = host->gap_block_list[block_id].size();
                    }
                }
                // else pos still in same consecutive gap
            }
            else
            {
                if (block_id && ((pos - gap_acc) % block_size) == block_size-1)  // switch blocks
                {
                    --block_id;
                    // is there are tailing gap in preceeding block overlapping with new pos?
                    if (host->gap_block_list[block_id].size() && host->gap_block_list[block_id].back().first + gap_acc == pos)
                    {
                            is_at_gap = true;
                            gap_id = host->gap_block_list[block_id].size() - 1;
                    }
                    else // no succeeding gap in same block, gap_id points to end of block
                        gap_id = host->gap_block_list[block_id].size();
                }
                else // no block switch
                {
                    if (gap_id && host->gap_block_list[block_id].size() && host->gap_block_list[block_id][gap_id-1].first + gap_acc - 1 == pos)
                    {
                        --gap_id;
                        is_at_gap = true;
                    }
                }
            }
            ungapped_view_pos = pos - gap_acc;
            return *this;
        }

        //!\brief Post-increment, returns previous iterator state (delegates to pre-increment).
        constexpr gap_decorator_anchor_block_iterator operator++(int) noexcept
        {
            gap_decorator_anchor_block_iterator cpy{*this};
            ++(*this);
            return cpy;
        }

        //!\brief Post-decrement, returns previous iterator state (delegates to pre-decrement).
        constexpr gap_decorator_anchor_block_iterator operator--(int) noexcept
        {
            gap_decorator_anchor_block_iterator cpy{*this};
            --(*this);
            return cpy;
        }
        //!\}

        /*!\name Reference/Dereference operators
         * \{
        */
        //!\brief Dereference operator returns a copy of the element currently pointed at.
        constexpr reference operator*() const noexcept
        {
            return (is_at_gap) ? static_cast<reference>(gap{})
                               : static_cast<reference>(host->ungapped_view[ungapped_view_pos]);
        }
        //!\}

        /*!\name Comparison operators
         * \brief Compares iterators by virtual position.
         * \{
         */
        constexpr friend bool operator==(gap_decorator_anchor_block_iterator const & lhs,
                                         gap_decorator_anchor_block_iterator const & rhs)
        {
            return lhs.pos == rhs.pos;
        }

        constexpr friend bool operator!=(gap_decorator_anchor_block_iterator const & lhs,
                                         gap_decorator_anchor_block_iterator const & rhs)
        {
            return lhs.pos != rhs.pos;
        }

        constexpr friend bool operator<(gap_decorator_anchor_block_iterator const & lhs,
                                        gap_decorator_anchor_block_iterator const & rhs)
        {
            return lhs.pos < rhs.pos;
        }

        constexpr friend bool operator>(gap_decorator_anchor_block_iterator const & lhs,
                                        gap_decorator_anchor_block_iterator const & rhs)
        {
            return lhs.pos > rhs.pos;
        }

        constexpr friend bool operator<=(gap_decorator_anchor_block_iterator const & lhs,
                                         gap_decorator_anchor_block_iterator const & rhs)
        {
            return lhs.pos <= rhs.pos;
        }

        constexpr friend bool operator>=(gap_decorator_anchor_block_iterator const & lhs,
                                         gap_decorator_anchor_block_iterator const & rhs)
        {
            return lhs.pos >= rhs.pos;
        }
        //!\}
    };

    //!\brief Short datatype for storing block and gap indices
    //using index_type = unsigned short int;

public:
    /*!\name Range-associated member types
     * \{
     */
    //!\brief The union type of the alphabet type and gap symbol type (see seqan3::gapped).
    using value_type = gapped<value_type_t<inner_type>>;
    //!\brief Use the value type as reference type because the underlying sequence must not be modified.
    using reference = value_type;
    //!\brief const_reference type equals reference type equals value type because the underlying sequence must not
    //!       be modified.
    using const_reference = reference;
    //!\brief The size_type of the underlying sequence.
    using size_type = size_type_t<inner_type>;
    //!\brief The difference type of the underlying sequence.
    using difference_type = difference_type_t<inner_type>;
    //!\brief The iterator type of this container (a bidirectional iterator).
    using iterator = gap_decorator_anchor_block_iterator;
    //!\brief The const_iterator equals the iterator type. Since no references are ever returned and thus the underlying
    //!        sequence cannot be modified through the iterator there is no need for const.
    using const_iterator = iterator;
    //!\}

    //!\brief The underlying ungapped range type.
    using unaligned_seq_type = inner_type;

    /*!\name Constructors, destructor and assignment.
     * \{
     */
    //!\brief Default constructor. Attention: all operations on a solely default constructed decorator,
    //!       except assigning a new range, are UB.
    constexpr gap_decorator_anchor_block() = default;
    //!\brief Copy constructor.
    constexpr gap_decorator_anchor_block(gap_decorator_anchor_block const &) = default;
    //!\brief Copy construction via assignment.
    constexpr gap_decorator_anchor_block & operator=(gap_decorator_anchor_block const &) = default;
    //!\brief Move constructor.
    constexpr gap_decorator_anchor_block(gap_decorator_anchor_block && rhs) = default;
    //!\brief Move assignment.
    constexpr gap_decorator_anchor_block & operator=(gap_decorator_anchor_block && rhs) = default;
    //!\brief Use default deconstructor.
    ~gap_decorator_anchor_block() = default;

    //!\brief Construct with the ungapped range type.
    template <typename other_range_t>
         requires !std::Same<other_range_t, gap_decorator_anchor_block> &&
                  std::Same<remove_cvref_t<other_range_t>, remove_cvref_t<inner_type>> &&
                  std::ranges::ViewableRange<other_range_t> // at end, otherwise it competes with the move ctor
    gap_decorator_anchor_block(other_range_t && range) : ungapped_view{std::view::all(std::forward<inner_type>(range))}
    {
        gap_block_list = gap_block_list_t(ungapped_view.size()/block_size + 1, gap_block_type{});
        gap_sums = gap_sums_type(ungapped_view.size()/block_size + 1, 0);
    } // TODO (@smehringer) only works for copyable views. Has to be changed once views are not required to be copyable anymore.
    // !\}

    /*!\brief Returns the total length of the aligned sequence.
     * \returns The total length of the aligned sequence (gaps included).
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    size_type size() const noexcept
    {
        if (gap_sums.size())
            return gap_sums.back() + ungapped_view.size();
        return ungapped_view.size();
    }

    /*!\name Aligned sequence modifications
     * \{
     */
    /*!\brief Insert a gap of length count at the aligned sequence iterator position.
     * \param it     Iterator indicating the gap start position in the aligned sequence.
     * \param count  Number of gap symbols to be inserted.
     * \returns      An iterator pointing to the start position of the insertion.
     *
     * ### Complexity
     *
     * Average and worst case (insertion before last gap): \f$O(k)\f$,
     * Best case (back insertion): \f$O(\log k)\f$.
     */
    iterator insert_gap(iterator const it, size_type const count = 1)
    {
        assert(ungapped_view.size());

        if (!count) // [[unlikely]]
            return it;

        size_type const pos = std::distance(begin(), it);
        assert(pos <= size());
        iterator it_tmp{*this, (pos) ? pos-1 : 0};

        size_type block_id = it.block_id;

        if (it.is_at_gap) // extension of current gap
            gap_block_list[it.block_id][it.gap_id].second += count;
        else if (it != begin() && it_tmp.is_at_gap) // extension of preceeding gap
        {
            gap_block_list[it_tmp.block_id][it_tmp.gap_id].second += count;
            block_id = it_tmp.block_id;
        }
        else // create new gap
            gap_block_list[it.block_id].insert(gap_block_list[it.block_id].begin() + it.gap_id, gap_t{pos - it.gap_acc, count});

        update_block_sums(block_id, count);

        return iterator{*this, pos};
    }

   /*!\brief Erase one gap symbol at the indicated iterator postion.
    * \param it     Iterator indicating the gap to be erased.
    * \returns      Iterator following the last removed element.
    * \throws seqan3::gap_erase_failure if character is no seqan3::gap.
    *
    * \details
    *
    * ### Complexity
    *
    * \f$O(\log k)\f$
    */
    iterator erase_gap(iterator const it)
    {
        assert(ungapped_view.size());

        // check if [it, it+gap_len[ covers [first, last[
        if ((*it) != gap{}) // [[unlikely]]
            throw gap_erase_failure("The range to be erased does not correspond to a consecutive gap.");

        auto end_it = std::next(it);
        return erase_gap(it, end_it);
    }

    /*!\brief Erase gap symbols at the iterator postions [first, last[.
     * \param[in]   first    The iterator pointing to the position where to start inserting gaps.
     * \param[in]   last     The iterator pointing to the position where to stop erasing gaps.
     * \returns     Iterator following the last removed element.
     * \throws seqan3::gap_erase_failure if [\p first, \p last[ does not correspond
     * to a consecutive range of seqan3::gap 's.
     *
     * \details
     *
     * ### Complexity
     *
     * \f$O(\log k)\f$
     */
    iterator erase_gap(iterator const first, iterator const last)
    {
        assert(ungapped_view.size());

        size_type const pos1 = std::distance(begin(), first);
        size_type const pos2 = std::distance(begin(), last);

        size_type block_id = first.block_id;
        size_type gap_id = first.gap_id;

        if (!first.is_at_gap || pos2 > gap_block_list[block_id][gap_id].first + first.gap_acc)
            throw gap_erase_failure("The range to be erased does not correspond to a consecutive gap.");

        // delete complete gap
        if (pos2 == size() || *last != gap{})
            gap_block_list[block_id].erase(gap_block_list[block_id].begin() + gap_id);

        else // decrease gap length
            gap_block_list[block_id][gap_id].second -= pos2 - pos1;

        update_block_sums(block_id, pos2-pos1, false);
        return iterator{*this, pos1};;
    }

    /*!\brief Assigns a new sequence of type seqan3:;gap_decorator_anchor_block::unaligned_seq_type to the decorator.
     * \param[in,out] dec       The decorator to modify.
     * \param[in]     unaligned The unaligned sequence to assign.
     */
    template <typename unaligned_seq_t> // generic template to use forwarding reference
    //!\cond
        requires std::Constructible<gap_decorator_anchor_block, unaligned_seq_t>
    //!\endcond
    friend void assign_unaligned(gap_decorator_anchor_block & dec, unaligned_seq_t && unaligned)
    {
        dec = unaligned;
    }
    //!\}

    /*!\name Iterators
     * \{
     */
    /*!\brief Returns an iterator to the first element of the container.
     * \returns Iterator to the first element.
     *
     * If the container is empty, the returned iterator will be equal to end().
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    iterator begin() const noexcept
    {
        return iterator{*this};
    }

    //!\copydoc begin()
    const_iterator cbegin() const noexcept
    {
        return const_iterator{*this};
    }

    /*!\brief Returns an iterator to the element following the last element of the decorator.
     * \returns Iterator to the behind last element.
     *
     * \attention This element acts as a placeholder; attempting to dereference it results in undefined behaviour.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    iterator end() const noexcept
    {
        return iterator{*this, size()};
    }

    //!\copydoc end()
    const_iterator cend() const noexcept
    {
        return const_iterator{*this, size()};
    }
    //!\}

    /*!\name Element access
     * \{
     */
    /*!\brief Return the i-th element as a reference.
     * \param i     The element to retrieve.
     * \returns     A reference of the gapped alphabet type.
     *
     * ### Complexity
     *
     * \f$O(\log k)\f$ where \f$k\f$ is the number of gaps.
     *
     * ### Exceptions
     *
     * Throws std::out_of_range exception if \p i is out of range.
     */
    reference at(size_type const i)
    {
        if (i >= size()) // [[unlikely]]
            throw std::out_of_range{"Trying to access element behind the last in gap_decorator."};
        return (*this)[i];
    }

    //!\copydoc at()
    const_reference at(size_type const i) const
    {
        if (i >= size()) // [[unlikely]]
            throw std::out_of_range{"Trying to access element behind the last in gap_decorator."};
        return (*this)[i];
    }

    /*!\brief Return the i-th element as a reference.
     * \param i     The element to retrieve.
     * \returns     A reference of the gapped alphabet type.
     *
     * This function delegates to an iterator seqan3::gap_decorator_anchor_block.
     *
     * ### Complexity
     *
     * \f$O(\log k)\f$ where \f$k\f$ is the number of gaps.
     *
     */
    constexpr reference operator[](size_type const i) const noexcept
    {
        return *iterator{*this, i};
    }
    //!\}

    /*!\name Comparison operators
     * \{
     */
    /*!\brief Compares two seqan3::gap_decorator_anchor_block 's by underlying sequence and gaps.
     * \param[in] lhs The left-hand side gap decorator to compare.
     * \param[in] rhs The right-hand side gap decorator to compare.
     * \returns A boolean flag indicating (in)equality of the aligned sequences.
     *
     * ### Complexity
     * Worst case: \f$O(n*\log k)\f$
     * Constant in case the decorators have not the same number of (consecutive) gaps.
     *
     * ### Exceptions
     *
     * No-throw guarantee. Does not modify the aligned sequences.
     */
    friend bool operator==(gap_decorator_anchor_block const & lhs, gap_decorator_anchor_block const & rhs) noexcept
    {
        if (lhs.size()  != rhs.size() || !std::ranges::equal(lhs.ungapped_view, rhs.ungapped_view))
            return false;

        std::vector<bool> gap_block_sim;
        std::function<bool(const gap_block_type &, const gap_block_type &)> block_comp =
            [](const gap_block_type & block_lhs, const gap_block_type & block_rhs)
            {return std::equal(block_lhs.begin(), block_lhs.end(), block_rhs.begin());};

        std::transform(lhs.gap_block_list.begin(), lhs.gap_block_list.end(), rhs.gap_block_list.begin(),
            std::back_inserter(gap_block_sim), block_comp);
        return (std::find(gap_block_sim.begin(), gap_block_sim.end(), false) == gap_block_sim.end()) ? true : false;
    }

    //!\copydoc operator==
    friend bool operator!=(gap_decorator_anchor_block const & lhs, gap_decorator_anchor_block const & rhs) noexcept
    {
        return !(lhs == rhs);
    }

    friend bool operator<(gap_decorator_anchor_block const & lhs, gap_decorator_anchor_block const & rhs) noexcept
    {
        auto lit = lhs.begin();
        auto rit = rhs.begin();

        while (lit != lhs.end() && rit != rhs.end() && *lit == *rit)
            ++lit, ++rit;

        if (rit == rhs.end())
            return false;           //  lhs == rhs, or rhs prefix of lhs
        else if (lit == lhs.end())
            return true;            // lhs prefix of rhs

        return *lit < *rit;
    }

    friend bool operator<=(gap_decorator_anchor_block const & lhs, gap_decorator_anchor_block const & rhs) noexcept
    {
        auto lit = lhs.begin();
        auto rit = rhs.begin();

        while (lit != lhs.end() && rit != rhs.end() && *lit == *rit)
            ++lit, ++rit;

        if (lit == lhs.end())
            return true;            // lhs == rhs, or lhs prefix of rhs
        else if (rit == rhs.end())
            return false;           // rhs prefix of lhs

        return *lit < *rit;
    }

    friend bool operator>(gap_decorator_anchor_block const & lhs, gap_decorator_anchor_block const & rhs) noexcept
    {
        return !(lhs <= rhs);
    }

    friend bool operator>=(gap_decorator_anchor_block const & lhs, gap_decorator_anchor_block const & rhs) noexcept
    {
        return !(lhs < rhs);
    }
    //!\}

private:
    //!\brief The gap type as a tuple storing position and length of a gap.
    using gap_t             = std::pair<size_type, size_type>;
    using gap_sums_type     = std::vector<size_type>;
    using gap_block_type    = std::vector<gap_t>;
    using gap_block_list_t  = std::vector<gap_block_type>;

    //!\brief Stores a (copy of a) view to the ungapped, underlying sequence.
    decltype(std::view::all(std::declval<inner_type &&>())) ungapped_view{};

    // block-wise accumulation of gap lengths, size doesn't change after sequence assignment
    std::vector<size_t> gap_sums{};
    // nested vector to store gaps by block, i.e. gap_block_list[block id] = gap_block_list
    // an eventually tailing gap will be stored in the last block for having consistent gap read behaviour in all blocks
    gap_block_list_t gap_block_list{};

    // Update block statistics
    void update_block_sums(size_type block_id, size_type count, bool add = true)
    {

        for (size_type id = block_id; id < gap_sums.size(); ++id)
            gap_sums[id] = (add) ? gap_sums[id] + count : gap_sums[id] - count;
    }
};

/*!\name Type deduction guides
 * \{
 */
//!\brief Ranges (not views!) always deduce to `const & range_type` since they are access-only anyway.
template <std::ranges::ViewableRange urng_t>
//!\cond
    requires !std::ranges::View<std::remove_reference_t<urng_t>>
//!\endcond
gap_decorator_anchor_block(urng_t && range) -> gap_decorator_anchor_block<std::remove_reference_t<urng_t> const &>;

//!\brief Views always deduce to their respective type because they are copied.
template <std::ranges::View urng_t>
gap_decorator_anchor_block(urng_t range) -> gap_decorator_anchor_block<urng_t>;
//!\}

} // namespace seqan

namespace seqan3::detail
{

//!\brief Type trait that declares any seqan3::gap_decorator_anchor_block to be **NOT a view**.
template <typename type>
constexpr int enable_view<seqan3::gap_decorator_anchor_block<type>> = 0;

template <typename type>
constexpr int enable_view<seqan3::gap_decorator_anchor_block<type> const> = 0;

} // namespace seqan3::detail

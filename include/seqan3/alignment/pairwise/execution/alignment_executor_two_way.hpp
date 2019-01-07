// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
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

/*!\file
 * \brief Provides seqan3::detail::alignment_executor_two_way.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <functional>
#include <optional>
#include <type_traits>

#include <seqan3/alignment/pairwise/execution/alignment_range.hpp>
#include <seqan3/alignment/pairwise/execution/execution_handler_sequential.hpp>
#include <seqan3/core/metafunction/range.hpp>
#include <seqan3/range/view/single_pass_input.hpp>
#include <seqan3/range/shortcuts.hpp>
#include <seqan3/std/ranges>

namespace seqan3::detail
{

/*!\brief A two way executor for pairwise alignments.
 * \ingroup execution
 * \tparam resource_t            The underlying range of sequence pairs to be computed; must model
 *                               std::ranges::ViewableRange and std::ranges::InputRange.
 * \tparam alignment_algorithm_t The alignment algorithm to be invoked on each sequence pair.
 * \tparam execution_handler_t   The execution handler managing the execution of the alignments.
 *
 * \details
 *
 * This alignment executor provides an additional buffer over the computed alignments to allow
 * a two-way execution flow. The alignment results can then be accessed in an order preserving manner using the
 * alignment_executor_two_way::bump() member function.
 */
template <std::ranges::ViewableRange resource_t,
          typename alignment_algorithm_t,
          typename execution_handler_t = execution_handler_sequential>
//!\cond
    requires std::ranges::InputRange<std::remove_reference_t<resource_t>>
//!\endcond
class alignment_executor_two_way
{
private:
    /*!\name Resource types
     * \{
     */
    //!\brief The underlying resource type.
    using resource_type       = detail::single_pass_input_view<resource_t>;
    //!\brief The value type of the resource.
    using resource_value_type = value_type_t<resource_type>;
    //!\}

    /*!\name Buffer types
     * \{
     */
    //!\brief The result of invoking the alignment instance.
    using buffer_value_type = typename alignment_algorithm_t::result_type;
    //!\brief The internal buffer.
    using buffer_type       = std::vector<buffer_value_type>;
    //!\brief The pointer type of the buffer.
    using buffer_pointer    = std::ranges::iterator_t<buffer_type>;
    //!\}

public:

    /*!\name Member types
     * \{
     */

    //!\brief The result type of invoking the alignment instance.
    using value_type      = buffer_value_type;
    //!\brief A reference to the alignment result.
    using reference       = std::add_lvalue_reference_t<value_type>;
    //!\brief The difference type for the buffer.
    using difference_type = typename buffer_type::difference_type;
    //!\}

    /*!\name Constructors, destructor and assignment
     * \brief The class is not copy-constructible or copy-assignable but allows move construction and assignment.
     * \{
     */
    alignment_executor_two_way()
    {
        init_buffer();
    }

    alignment_executor_two_way(alignment_executor_two_way const &) = delete;
    alignment_executor_two_way(alignment_executor_two_way &&) = default;
    alignment_executor_two_way & operator=(alignment_executor_two_way const &) = delete;
    alignment_executor_two_way & operator=(alignment_executor_two_way && ) = default;
    ~alignment_executor_two_way() = default;

    //!\brief Constructs this executor with the passed range of alignment instances.
    alignment_executor_two_way(resource_t && resrc,
                               alignment_algorithm_t fn) :
        resource{view::single_pass_input(std::forward<resource_t>(resrc))},
        kernel{fn}
    {
        init_buffer();
    }
    //!}

    /*!\name Modifiers
     * \{
     */
    /*!\brief Resets the underlying resource.
     * \param[in] resrc The new range of sequence pairs to compute alignments for.
     *
     * \attention Setting a new resource will discard the currently buffered alignment results immediately.
     */
    void reset_resource(resource_t && resrc)
    {
        resource = std::forward<resource_t>(resrc) | view::single_pass_input;
        setg(seqan3::end(buffer), seqan3::end(buffer));
    }

    /*!\brief Resets the underlying algorithm to compute the alignments.
     * \param[in] fn The new algorithm to compute the alignments with.
     *
     * \attention Setting a new alignment function will discard the currently buffered alignment results immediately.
     */
    void reset_algorithm(alignment_algorithm_t fn)
    {
        kernel = fn;
        setg(seqan3::end(buffer), seqan3::end(buffer));
    }
    //!\}

    /*!\name Get area
     * \{
     */
    /*!\brief Returns the current alignment result in the buffer and advances the buffer to the next position.
     * \returns A std::optional that either contains a reference to the underlying value or is empty, iff the
     *          underlying resource has been completely consumed.
     *
     * \details
     *
     * If there is no available input in the result buffer anymore, this function triggers an underflow to fill
     * the buffer with the next alignments.
     *
     * ### Exception
     *
     * Throws std::bad_function_call if the algorithm was not set.
     */
    std::optional<std::reference_wrapper<value_type>> bump()
    {
        if (gptr == buffer_pointer{} || in_avail() == 0)
        {
            if (underflow() == eof)
            {
                return {};
            }
        }
        return {std::ref(*gptr++)};
    }

    //!\brief Returns the remaining number of elements in the buffer, that are not read yet.
    constexpr size_t in_avail() const noexcept
    {
        return egptr - gptr;
    }
    //!\}

    /*!\name Miscellaneous
     * \{
     */
    //!\brief Checks whether the end of the input resource was reached.
    bool is_eof() noexcept
    {
        return seqan3::begin(resource) == seqan3::end(resource);
    }
    //!\}

private:

    /*!\name Get area
     * \{
     */

    //\brief Sets the buffer pointer.
    void setg(buffer_pointer beg, buffer_pointer end)
    {
        gptr = beg;
        egptr = end;
    }

    //!\brief Refills the buffer with new alignment results.
    size_t underflow()
    {
        if (gptr < egptr)  // Case: buffer not completely consumed
            return in_avail();

        if (is_eof())  // Case: reached end of resource.
            return eof;

        // Reset the get pointer.
        setg(seqan3::begin(buffer), seqan3::end(buffer));

        // Apply the alignment execution.
        // TODO: Adapt for async behavior for parallel execution handler.
        size_t count = 0;
        for (auto resource_iter = seqan3::begin(resource);
             count < in_avail() && !is_eof(); ++count, ++resource_iter, ++gptr)
        {
            auto const & [first_seq, second_seq] = *resource_iter;
            exec_handler.execute(kernel, first_seq, second_seq, [this](auto && res){ *gptr = std::move(res); });
        }

        // Update the available get position if the buffer was consumed completely.
        setg(seqan3::begin(buffer), seqan3::begin(buffer) + count);

        return in_avail();
    }
    //!\}

    /*!\name Miscellaneous
     * \{
     */

    //!\brief Initializes the underlying buffer.
    void init_buffer()
    {
        buffer.resize(1);
        setg(seqan3::end(buffer), seqan3::end(buffer));
    }
    //!\}

private:

    //!\brief Indicates the end-of-stream.
    static constexpr size_t eof{-1};

    //!\brief The execution policy.
    execution_handler_t exec_handler{};

    //!\brief The underlying resource containing the alignment instances.
    resource_type resource{};  // view or lvalue ref
    //!\brief Selects the correct alignment to execute.
    alignment_algorithm_t kernel{};

    //!\brief The buffer storing the alignment results.
    buffer_type buffer{};
    //!\brief The get pointer in the buffer.
    buffer_pointer gptr{};
    //!\brief The end get pointer in the buffer.
    buffer_pointer egptr{};
};

/*!\name Type deduction guides
 * \relates seqan3::detail::alignment_executor_two_way
 * \{
 */
template <typename resource_rng_t, typename func_t>
alignment_executor_two_way(resource_rng_t &&, func_t) ->
    alignment_executor_two_way<resource_rng_t, func_t, execution_handler_sequential>;

//!\}
} // namespace seqan3::detail

// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::alignment_executor_two_way.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <functional>
#include <optional>
#include <type_traits>

#include <seqan3/alignment/pairwise/alignment_range.hpp>
#include <seqan3/alignment/pairwise/execution/execution_handler_parallel.hpp>
#include <seqan3/alignment/pairwise/execution/execution_handler_sequential.hpp>
#include <seqan3/core/parallel/execution.hpp>
#include <seqan3/core/type_traits/range.hpp>
#include <seqan3/range/shortcuts.hpp>
#include <seqan3/range/view/single_pass_input.hpp>
#include <seqan3/range/view/view_all.hpp>
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
 * a two-way execution flow. The alignment results can then be accessed in an order-preserving manner using the
 * alignment_executor_two_way::bump() member function.
 */
template <std::ranges::ViewableRange resource_t,
          typename alignment_algorithm_t,
          typename execution_handler_t = execution_handler_sequential>
//!\cond
    requires std::ranges::ForwardRange<resource_t> &&
             std::CopyConstructible<alignment_algorithm_t>
//!\endcond
class alignment_executor_two_way
{
private:
    /*!\name Resource types
     * \{
     */
    //!\brief The underlying resource type.
    using resource_type       = decltype(view::single_pass_input(std::view::zip(std::declval<resource_t>(),
                                                                                std::view::iota(0))));
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
    alignment_executor_two_way() = delete;                                               //!< This is a move-only type.
    alignment_executor_two_way(alignment_executor_two_way const &) = delete;             //!< This is a move-only type.
    alignment_executor_two_way(alignment_executor_two_way &&) = default;                 //!< Defaulted
    alignment_executor_two_way & operator=(alignment_executor_two_way const &) = delete; //!< This is a move-only type.
    alignment_executor_two_way & operator=(alignment_executor_two_way && ) = default;    //!< Defaulted
    ~alignment_executor_two_way() = default;                                             //!< Defaulted

    /*!\brief Constructs this executor with the passed range of alignment instances.
     * \param[in] resrc The underlying resource containing the sequence pairs to align.
     * \param[in] fn    The alignment kernel to invoke on the sequences pairs.
     *
     * \details
     *
     * Forwards the resource range as a zipped view with an index view to provide internal ids for the alignments.
     * If the execution handler is parallel it allocates a buffer of the size of the given resource range.
     * Otherwise the buffer size is 1.
     */
    alignment_executor_two_way(resource_t resrc,
                               alignment_algorithm_t fn) :
        resource{std::view::zip(std::forward<resource_t>(resrc), std::view::iota(0))},
        kernel{std::move(fn)}
    {
        if constexpr (std::Same<execution_handler_t, execution_handler_parallel>)
            init_buffer(std::ranges::distance(resrc));
        else
            init_buffer(1);
    }

    /*!\brief Constructs this executor with the passed range of alignment instances.
     * \param[in] resrc The underlying resource containing the sequence pairs to align.
     * \param[in] fn    The alignment kernel to invoke on the sequences pairs.
     * \param[in] exec  The execution policy to use.
     *
     * \details
     *
     * This constructor is used to type deduce the execution handler type from the passed execution policy.
     * The given execution policy is ignored.
     */
    template <typename exec_policy_t>
        requires is_execution_policy_v<exec_policy_t>
    alignment_executor_two_way(resource_t resrc,
                               alignment_algorithm_t fn,
                               exec_policy_t const & SEQAN3_DOXYGEN_ONLY(exec)) :
        alignment_executor_two_way{std::move(resrc), std::move(fn)}
    {
        static_assert(!std::Same<exec_policy_t, parallel_unsequenced_policy>,
                      "Parallel unsequenced execution not supported!");
        static_assert(!std::Same<exec_policy_t, unsequenced_policy>,
                      "Unsequenced execution not supported!");
    }
    //!}

    /*!\name Get area
     * \{
     */
    /*!\brief Returns the current alignment result in the buffer and advances the buffer to the next position.
     * \returns A std::optional that either contains a reference to the underlying value or is empty iff the
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
    std::optional<value_type> bump()
    {
        if (gptr == buffer_pointer{} || in_avail() == 0)
        {
            if (underflow() == eof)
            {
                return {std::nullopt};
            }
        }
        return {std::move(*gptr++)};
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
        return std::ranges::begin(resource) == std::ranges::end(resource);
    }
    //!\}

private:

    /*!\name Get area
     * \{
     */

    //!\brief Sets the buffer pointer.
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
        setg(std::ranges::begin(buffer), std::ranges::end(buffer));

        // Apply the alignment execution.
        size_t count = 0;
        size_t buffer_limit = in_avail();
        for (auto resource_iter = std::ranges::begin(resource);
             count < buffer_limit && !is_eof(); ++count, ++resource_iter, ++gptr)
        {
            auto && [tpl, idx] = *resource_iter;
            auto && [first_seq, second_seq] = tpl;
            buffer_pointer write_to = gptr;
            exec_handler.execute(kernel,
                                 idx,
                                 first_seq | view::all,
                                 second_seq | view::all,
                                 [write_to] (auto && res) { *write_to = std::move(res); });
        }

        exec_handler.wait();

        // Update the available get position if the buffer was consumed completely.
        setg(std::ranges::begin(buffer), std::ranges::begin(buffer) + count);

        return in_avail();
    }
    //!\}

    /*!\name Miscellaneous
     * \{
     */

    /*!\brief Initialises the underlying buffer.
     * \param size The initial size of the buffer.
     */
    void init_buffer(size_t const size)
    {
        buffer.resize(size);
        setg(std::ranges::end(buffer), std::ranges::end(buffer));
    }
    //!\}

    //!\brief Indicates the end-of-stream.
    static constexpr size_t eof{std::numeric_limits<size_t>::max()};

    //!\brief The execution policy.
    execution_handler_t exec_handler{};

    //!\brief The underlying resource containing the alignment instances.
    resource_type resource{};
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

//!\brief Deduce the type from the provided arguments and set the sequential execution handler.
template <typename resource_rng_t, typename func_t>
alignment_executor_two_way(resource_rng_t &&, func_t) ->
    alignment_executor_two_way<resource_rng_t, func_t, execution_handler_sequential>;

//!\brief Deduce the type from the provided arguments and set the sequential execution handler.
template <typename resource_rng_t, typename func_t, typename exec_policy_t>
    requires is_execution_policy_v<exec_policy_t>
alignment_executor_two_way(resource_rng_t &&, func_t, exec_policy_t const &) ->
    alignment_executor_two_way<resource_rng_t,
                               func_t,
                               std::conditional_t<std::Same<exec_policy_t, parallel_policy>,
                                                  execution_handler_parallel,
                                                  execution_handler_sequential>>;

//!\}
} // namespace seqan3::detail

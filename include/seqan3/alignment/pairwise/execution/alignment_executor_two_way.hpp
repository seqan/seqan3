// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
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
#include <seqan3/alignment/pairwise/detail/type_traits.hpp>
#include <seqan3/alignment/pairwise/execution/execution_handler_parallel.hpp>
#include <seqan3/alignment/pairwise/execution/execution_handler_sequential.hpp>
#include <seqan3/core/parallel/execution.hpp>
#include <seqan3/core/type_traits/range.hpp>
#include <seqan3/range/views/chunk.hpp>
#include <seqan3/range/views/type_reduce.hpp>
#include <seqan3/range/views/zip.hpp>
#include <seqan3/std/ranges>

namespace seqan3::detail
{

/*!\brief A two way executor for pairwise alignments.
 * \ingroup execution
 * \tparam resource_t            The underlying range of sequence pairs to be computed; must model
 *                               std::ranges::viewable_range and std::ranges::input_range.
 * \tparam alignment_algorithm_t The alignment algorithm to be invoked on each sequence pair.
 * \tparam execution_handler_t   The execution handler managing the execution of the alignments.
 *
 * \details
 *
 * This alignment executor provides an additional buffer over the computed alignments to allow
 * a two-way execution flow. The alignment results can then be accessed in an order-preserving manner using the
 * alignment_executor_two_way::bump() member function.
 */
template <std::ranges::viewable_range resource_t,
          typename alignment_algorithm_t,
          typename execution_handler_t = execution_handler_sequential>
//!\cond
    requires std::ranges::forward_range<resource_t> &&
             std::copy_constructible<alignment_algorithm_t>
//!\endcond
class alignment_executor_two_way
{
private:
    /*!\name Resource types
     * \{
     */
    //!\brief The underlying resource type augmented with an index and a chunked view.
    using chunked_resource_type = typename chunked_indexed_sequence_pairs<resource_t>::type;
    //!\brief The iterator over the underlying resource.
    using chunked_resource_iterator = std::ranges::iterator_t<chunked_resource_type>;
    //!\}

    /*!\name Buffer types
     * \{
     */
    //!\brief The result of invoking the alignment instance.
    using alignment_result_type = typename alignment_algorithm_t::result_type;
    //!\brief The actual alignment result.
    using buffer_value_type = std::ranges::range_value_t<alignment_result_type>;
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
     * \brief The class is move-only, i.e. it is not copy-constructible or copy-assignable.
     * \{
     */
    //!\brief Deleted default constructor because this class manages an external resource.
    alignment_executor_two_way() = delete;
    //!\brief This class provides unique ownership over the managed resource and is therefor not copyable.
    alignment_executor_two_way(alignment_executor_two_way const &) = delete;

    /*!\brief Move constructs the resource of the other executor.
     * \param[in] other The other alignment executor (prvalue) to move from.
     *
     * \details
     *
     * Handling the move of the underlying resource, respectively result buffer, requires some non-default operations.
     * The iterator holding the current state of the executor must be reinitailised after the resource and buffer have
     * been moved.
     *
     * ### Exception
     *
     * no-throw guarantee.
     *
     * ### Complexity
     *
     * Constant if the underlying resource type models std::ranges::random_access_range, otherwise linear.
     */
    alignment_executor_two_way(alignment_executor_two_way && other) noexcept
    {
        move_initialise(std::move(other));
    }

    //!\brief This class provides unique ownership over the managed resource and is therefor not copyable.
    alignment_executor_two_way & operator=(alignment_executor_two_way const &) = delete;

    //!\brief Move assigns from the resource of another executor.
    //!\copydetails seqan3::detail::alignment_executor_two_way::alignment_executor_two_way(alignment_executor_two_way && other)
    alignment_executor_two_way & operator=(alignment_executor_two_way && other)
    {
        move_initialise(std::move(other));
        return *this;
    }

    //!\brief Defaulted.
    ~alignment_executor_two_way() = default;

    /*!\brief Constructs this executor with the passed range of alignment instances.
     * \tparam exec_policy_t The type of the execution policy; seqan3::is_execution_policy must return `true`. Defaults
     *                       to seqan3::sequenced_policy.
     * \param[in] resrc The underlying resource containing the sequence pairs to align.
     * \param[in] fn The alignment kernel to invoke on the sequences pairs.
     * \param[in] chunk_size The number of sequence pairs to invoke the alignment algorithm with.
     * \param[in] exec Optional execution policy to use. Defaults to seqan3::seq.
     *
     * \throws std::invalid_argument if the chunk size is less than 1.
     *
     * \details
     *
     * Forwards the resource range as a zipped view with an index view to provide internal ids for the alignments.
     * If the execution handler is parallel, it allocates a buffer of the size of the given resource range.
     * Otherwise the buffer size is 1.
     */
    template <typename exec_policy_t = sequenced_policy>
    //!\cond
        requires is_execution_policy_v<exec_policy_t>
    //!\endcond
    alignment_executor_two_way(resource_t resrc,
                               alignment_algorithm_t fn,
                               size_t chunk_size = 1u,
                               exec_policy_t const & SEQAN3_DOXYGEN_ONLY(exec) = seq) :
        kernel{std::move(fn)},
        _chunk_size{chunk_size}
    {

        if (chunk_size == 0u)
            throw std::invalid_argument{"The chunk size must be greater than 0."};

        chunked_resource = views::zip(std::forward<resource_t>(resrc), std::views::iota(0)) | views::chunk(_chunk_size);
        chunked_resource_it = chunked_resource.begin();

        if constexpr (std::same_as<execution_handler_t, execution_handler_parallel>)
            init_buffer(std::ranges::distance(resrc));
        else
            init_buffer(_chunk_size);
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
        if (underflow() == eof)
            return {std::nullopt};

        assert(gptr < egptr);
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
        return chunked_resource_it == std::ranges::end(chunked_resource);
    }

    //!\brief Returns the selected chunk size.
    constexpr size_t chunk_size() const noexcept
    {
        return _chunk_size;
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

        using std::get;

        // Apply the alignment execution.
        size_t count = 0;
        size_t buffer_limit = in_avail();

        for (size_t current_chunk_size = std::ranges::distance(*chunked_resource_it);
             count < buffer_limit && !is_eof();
             count += current_chunk_size, gptr += current_chunk_size, ++chunked_resource_it)
        {
            auto && current_chunk = std::ranges::iter_move(chunked_resource_it);
            current_chunk_size = std::ranges::distance(current_chunk);
            assert(in_avail() >= current_chunk_size);

            exec_handler.execute(kernel, std::move(current_chunk), [write_to = gptr] (auto && alignment_results)
            {
                std::ranges::move(alignment_results, write_to);
            });
        }

        exec_handler.wait();

        // Update the available get positions to cover the part of the buffer that was filled.
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

    //!\brief Helper function to move initialise `this` from `other`.
    //!\copydetails seqan3::detail::alignment_executor_two_way::alignment_executor_two_way(alignment_executor_two_way && other)
    void move_initialise(alignment_executor_two_way && other) noexcept
    {
        kernel = std::move(other.kernel);
        _chunk_size = std::move(other._chunk_size);
        // Get the old resource position.
        std::ptrdiff_t old_resource_pos = std::ranges::distance(other.chunked_resource.begin(),
                                                                other.chunked_resource_it);
        // Move the resource and set the iterator state accordingly.
        chunked_resource = std::move(other.chunked_resource);
        chunked_resource_it = std::ranges::next(chunked_resource.begin(), old_resource_pos);

        // Get the old get pointer positions.
        std::ptrdiff_t old_gptr_pos = other.gptr - buffer.begin();
        std::ptrdiff_t old_egptr_pos = other.egptr - buffer.begin();
        // Move the buffer and set the get pointer accordingly.
        buffer = std::move(other.buffer);
        setg(buffer.begin() + old_gptr_pos, buffer.begin() + old_egptr_pos);
    }
    //!\}

    //!\brief Indicates the end-of-stream.
    static constexpr size_t eof{std::numeric_limits<size_t>::max()};

    //!\brief The execution policy.
    execution_handler_t exec_handler{};

    //!\brief The underlying resource containing the alignment instances.
    chunked_resource_type chunked_resource{};
    //!\brief The iterator over the resource that stores the current state of the executor.
    chunked_resource_iterator chunked_resource_it{};
    //!\brief Selects the correct alignment to execute.
    alignment_algorithm_t kernel{};

    //!\brief The buffer storing the alignment results.
    buffer_type buffer{};
    //!\brief The get pointer in the buffer.
    buffer_pointer gptr{};
    //!\brief The end get pointer in the buffer.
    buffer_pointer egptr{};
    //!\brief The size of the chunks to call the stored algorithm with.
    size_t _chunk_size{};
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
alignment_executor_two_way(resource_rng_t &&, func_t, size_t, exec_policy_t const &) ->
    alignment_executor_two_way<resource_rng_t,
                               func_t,
                               std::conditional_t<std::same_as<exec_policy_t, parallel_policy>,
                                                  execution_handler_parallel,
                                                  execution_handler_sequential>>;

//!\}
} // namespace seqan3::detail

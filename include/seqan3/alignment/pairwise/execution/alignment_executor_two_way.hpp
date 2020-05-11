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
 * \tparam value_t               The value type to buffer results generated from the algorithm; must model
 *                               std::semiregular.
 * \tparam execution_handler_t   The execution handler managing the execution of the alignments.
 *
 * \details
 *
 * This alignment executor provides an additional buffer over the computed alignments to allow
 * a two-way execution flow. The alignment results can then be accessed in an order-preserving manner using the
 * alignment_executor_two_way::bump() member function.
 *
 * ### Bucket structure
 *
 * Since it is not clear how many results a single invocation of the given algorithm produces the buffered results
 * are placed into buckets. The number of available buckets is determined by the execution policy. In sequential
 * execution mode only one bucket is available and only one invocation is buffered at a time. In the parallel execution,
 * a bucket is allocated for every element of the underlying resource.
 */
template <std::ranges::viewable_range resource_t,
          typename alignment_algorithm_t,
          std::semiregular value_t,
          typename execution_handler_t = execution_handler_sequential>
//!\cond
    requires std::ranges::forward_range<resource_t> && std::copy_constructible<alignment_algorithm_t>
//!\endcond
class alignment_executor_two_way
{
private:
    /*!\name Resource types
     * \{
     */
    //!\brief The underlying resource type.
    using resource_type = std::ranges::all_view<resource_t>;
    //!\brief The iterator over the underlying resource.
    using resource_iterator_type = std::ranges::iterator_t<resource_type>;
    //!\}

    /*!\name Buffer types
     * \{
     */
    //!\brief The internal buffer.
    using buffer_type = std::vector<value_t>;
    //!\brief The type of the bucket container.
    using bucket_type = std::vector<buffer_type>;
    //!\brief The iterator type of the buffer.
    using buffer_iterator_type = std::ranges::iterator_t<buffer_type>;
    //!\brief The iterator type over the bucket container.
    using bucket_iterator_type = std::ranges::iterator_t<bucket_type>;
    //!\}

    //!\brief Return status for seqan3::detail::alignment_executor_two_way::underflow.
    enum underflow_status
    {
        non_empty_buffer, //!< The buffer is not fully consumed yet and contains at least one element.
        empty_buffer, //!< The buffer is empty after calling underflow.
        end_of_resource //!< The end of the resource was reached.
    };

public:

    /*!\name Member types
     * \{
     */
    //!\brief The result type of invoking the alignment instance.
    using value_type      = value_t;
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
     *
     * \param[in] resource The underlying resource containing the sequence pairs to align.
     * \param[in] fn The alignment kernel to invoke on the sequences pairs.
     * \param[in] buffer_value A dummy object to deduce the type of the underlying buffer value.
     * \param[in] exec Optional execution policy to use. Defaults to seqan3::seq.
     *
     * \throws std::invalid_argument if the chunk size is less than 1.
     *
     * \details
     *
     * If the execution handler is parallel, it allocates a buffer of the size of the given resource range.
     * Otherwise the buffer size is 1.
     * Also note that the third argument is used for deducing the type of the underlying buffer value and is otherwise
     * not used in the context of the class' construction.
     */
    template <typename exec_policy_t = sequenced_policy>
    //!\cond
        requires is_execution_policy_v<exec_policy_t>
    //!\endcond
    alignment_executor_two_way(resource_t resource,
                               alignment_algorithm_t fn,
                               value_t const SEQAN3_DOXYGEN_ONLY(buffer_value) = value_t{},
                               exec_policy_t const & SEQAN3_DOXYGEN_ONLY(exec) = seq) :
        resource{std::views::all(resource)},
        resource_it{std::ranges::begin(this->resource)},
        kernel{std::move(fn)}
    {
        if constexpr (std::same_as<execution_handler_t, execution_handler_parallel>)
            bucket_vector_size = std::ranges::distance(resource);

        bucket_vector.resize(bucket_vector_size);
        bucket_iterator = bucket_vector.end();
        bucket_end = bucket_vector.end();
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
        underflow_status status;
        // Each invocation of the algorithm might produce zero results (e.g. a search might not find a query)
        // this repeats the algorithm until it produces the first result or the input resource was consumed.
        do { status = underflow(); } while (status == underflow_status::empty_buffer);

        if (status == underflow_status::end_of_resource)
            return {std::nullopt};

        assert(status == underflow_status::non_empty_buffer);
        assert(buffer_iterator != bucket_iterator->end());

        std::optional<value_type> element = std::ranges::iter_move(buffer_iterator);
        next_buffer_iterator(); // Go to next buffered value
        return element;
    }
    //!\}

    /*!\name Miscellaneous
     * \{
     */
    //!\brief Checks whether the end of the input resource was reached.
    bool is_eof() noexcept
    {
        return resource_it == std::ranges::end(resource);
    }
    //!\}

private:

    /*!\name Get area
     * \{
     */

    //!\brief Fills pre-assigned buckets (one bucket = results of one alignment) with new alignment results.
    underflow_status underflow()
    {
        if (!is_buffer_empty())  // Not everything consumed yet.
            return underflow_status::non_empty_buffer;

        if (is_eof())  // Case: reached end of resource.
            return underflow_status::end_of_resource;

        // Reset the buckets and the bucket iterator
        reset_buckets();

        // Execute the algorithm (possibly asynchronous) and fill the buckets in this pre-assigned order.
        for (bucket_end = bucket_iterator; bucket_end != bucket_vector.end() && !is_eof(); ++bucket_end, ++resource_it)
        {
            exec_handler.execute(kernel, *resource_it, [target_bucket_iterator = bucket_end] (auto && alignment_result)
            {
                target_bucket_iterator->push_back(std::move(alignment_result));
            });
        }

        exec_handler.wait();

        // Move the results iterator to the next available result. (This skips empty results of the algorithm)
        find_next_non_empty_bucket();

        if (is_buffer_empty())
            return underflow_status::empty_buffer;

        return underflow_status::non_empty_buffer;
    }
    //!\}

    /*!\name Miscellaneous
     * \{
     */

    /*!\brief Whether the internal buffer is empty.
     * \returns `true` if all elements of the internal buffer have been consumed, otherwise `false`.
     */
    bool is_buffer_empty() const
    {
        return bucket_iterator == bucket_end;
    }

    /*!\brief Resets the buckets.
     *
     * \details
     *
     * Clears all buckets and sets the bucket iterator to the first iterator, such that the allocated memory for each
     * bucket can be reused between invocations of seqan3::detail::alignment_executor_two_way::underflow.
     */
    void reset_buckets()
    {
        // Clear all buckets
        for (auto & bucket : bucket_vector)
            bucket.clear();

        // Reset the iterator over the buckets.
        bucket_iterator = bucket_vector.begin();
    }

    /*!\brief Finds the first non-empty bucket starting from the current bucket iterator.
     *
     * \details
     *
     * Finds the first non-empty bucket and sets the buffer iterator to the first element of this bucket.
     * If all buckets are empty, then the bucket iterator is set to the end of the buckets container and the buffer
     * iterator is not modified.
     */
    void find_next_non_empty_bucket()
    {
        assert(bucket_iterator <= bucket_end);
        // find first buffered bucket that contains at least one element
        bucket_iterator = std::find_if(bucket_iterator, bucket_end, [] (auto const & buffer)
        {
            return !buffer.empty();
        });

        if (bucket_iterator != bucket_end)
            buffer_iterator = bucket_iterator->begin();
    }

    /*!\brief Moves the buffer iterator to the next available element.
     *
     * \details
     *
     * If the current bucket is consumed, then the bucket iterator is incremented and the next non-empty bucket is found
     * by calling seqan3::detail::alignment_executor_two_way::find_next_non_empty_bucket.
     */
    void next_buffer_iterator()
    {
        if (++buffer_iterator == bucket_iterator->end())
        {
            ++bucket_iterator;
            find_next_non_empty_bucket();
        }
    }

    //!\brief Helper function to move initialise `this` from `other`.
    //!\copydetails seqan3::detail::alignment_executor_two_way::alignment_executor_two_way(alignment_executor_two_way && other)
    void move_initialise(alignment_executor_two_way && other) noexcept
    {
        kernel = std::move(other.kernel);
        bucket_vector_size = std::move(other.bucket_vector_size);
        // Get the old resource position.
        auto old_resource_position = std::ranges::distance(std::ranges::begin(other.resource),
                                                           other.resource_it);
        // Move the resource and set the iterator state accordingly.
        resource = std::move(other.resource);
        resource_it = std::ranges::next(std::ranges::begin(resource), old_resource_position);

        // Get the old get pointer positions.
        auto bucket_iterator_position = other.bucket_iterator - other.bucket_vector.begin();
        auto bucket_end_position = other.bucket_end - other.bucket_vector.begin();

        std::ptrdiff_t buffer_iterator_position = 0;
        if (bucket_iterator_position != bucket_end_position)
            buffer_iterator_position = other.buffer_iterator - other.bucket_iterator->begin();

        // Move the buffer and set the get pointer accordingly.
        bucket_vector = std::move(other.bucket_vector);
        bucket_iterator = bucket_vector.begin() + bucket_iterator_position;
        bucket_end = bucket_vector.begin() + bucket_end_position;

        if (bucket_iterator_position != bucket_end_position)
            buffer_iterator = bucket_iterator->begin() + buffer_iterator_position;
    }
    //!\}

    //!\brief The execution policy.
    execution_handler_t exec_handler{};

    //!\brief The underlying resource containing the alignment instances.
    resource_type resource{};
    //!\brief The iterator over the resource that stores the current state of the executor.
    resource_iterator_type resource_it{};
    //!\brief Selects the correct alignment to execute.
    alignment_algorithm_t kernel{};

    //!\brief The buffer storing the alignment results.
    bucket_type bucket_vector{};

    //!\brief The iterator pointing to the current bucket.
    bucket_iterator_type bucket_iterator{};
    //!\brief The iterator pointing behind the last bucket (must not be the end of the bucket_vector).
    bucket_iterator_type bucket_end{};

    //!\brief The buffer iterator pointing to the current result to be processed.
    buffer_iterator_type buffer_iterator{};
    //!\brief The end get pointer in the buffer.
    size_t bucket_vector_size{1};
};

/*!\name Type deduction guides
 * \relates seqan3::detail::alignment_executor_two_way
 * \{
 */

//!\brief Deduce the type from the provided arguments and set the sequential execution handler.
template <typename resource_rng_t, typename func_t, std::semiregular value_t>
alignment_executor_two_way(resource_rng_t &&, func_t, value_t const &) ->
    alignment_executor_two_way<resource_rng_t, func_t, value_t, execution_handler_sequential>;

//!\}
} // namespace seqan3::detail

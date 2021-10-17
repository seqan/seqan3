// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::algorithm_executor_blocking.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <functional>
#include <optional>
#include <seqan3/std/ranges>
#include <type_traits>

#include <seqan3/core/algorithm/detail/execution_handler_parallel.hpp>
#include <seqan3/core/algorithm/detail/execution_handler_sequential.hpp>

namespace seqan3::detail
{

/*!\brief A blocking algorithm executor for algorithms.
 * \ingroup core_algorithm
 * \tparam resource_t The underlying range of sequence pairs to be computed; must model std::ranges::viewable_range
 *                    and std::ranges::forward_range.
 * \tparam algorithm_t The algorithm to be invoked on the elements of the given resource; must model std::semiregular.
 * \tparam algorithm_result_t The result type generated by the algorithm; must model std::semiregular.
 * \tparam execution_handler_t The execution handler managing the execution of the algorithms.
 *
 * \details
 *
 * This blocking algorithm executor provides an additional buffer over the computed algorithm results to allow
 * a two-way execution flow. The results can then be accessed in an order-preserving manner using the
 * seqan3::detail::algorithm_executor_blocking::next_result() member function.
 *
 * ### Invocation
 *
 * The blocking algorithm executor invokes the algorithm on each element of the given resource. It passes as second
 * argument a callback function that stores the results of the algorithm in a pre-assigned buffer location.
 * The algorithm must therefor model std::invocable with the first argument type being std::ranges::range_reference_t
 * of the `resource_t` and the second argument type being convertible to std::function<void(algorithm_result_t)>.
 *
 * ### Result buffer
 *
 * Since it is not clear how many results a single invocation of the given algorithm produces the buffered results
 * are placed into buckets. The number of available buckets is determined by the execution policy. In sequential
 * execution mode only one bucket is available and only one invocation is buffered at a time. In the parallel execution,
 * a bucket is allocated for every element of the underlying resource.
 */
template <std::ranges::viewable_range resource_t,
          std::semiregular algorithm_t,
          std::semiregular algorithm_result_t,
          typename execution_handler_t = execution_handler_sequential>
//!\cond
    requires std::ranges::forward_range<resource_t> &&
             std::invocable<algorithm_t, std::ranges::range_reference_t<resource_t>,
                                         std::function<void(algorithm_result_t)>>
//!\endcond
class algorithm_executor_blocking
{
private:
    /*!\name Resource types
     * \{
     */
    //!\brief The underlying resource type.
    using resource_type = std::views::all_t<resource_t>;
    //!\brief The iterator over the underlying resource.
    using resource_iterator_type = std::ranges::iterator_t<resource_type>;
    //!\brief The difference type over the underlying resource.
    using resource_difference_type = std::iter_difference_t<resource_iterator_type>;
    //!\}

    /*!\name Buffer types
     * \{
     */
    //!\brief The type of a bucket storing the results produced by a single algorithm invocation.
    using bucket_type = std::vector<algorithm_result_t>;
    //!\brief The iterator type of a bucket.
    using bucket_iterator_type = std::ranges::iterator_t<bucket_type>;
    //!\brief The type of the buffer.
    using buffer_type = std::vector<bucket_type>;
    //!\brief The iterator type of the buffer.
    using buffer_iterator_type = std::ranges::iterator_t<buffer_type>;
    //!\}

    //!\brief Return status for seqan3::detail::algorithm_executor_blocking::fill_buffer.
    enum fill_status
    {
        non_empty_buffer, //!< The buffer is not fully consumed yet and contains at least one element.
        empty_buffer, //!< The buffer is empty after calling fill_buffer.
        end_of_resource //!< The end of the resource was reached.
    };

public:
    /*!\name Constructors, destructor and assignment
     * \brief The class is move-only, i.e. it is not copy-constructible or copy-assignable.
     * \{
     */
    //!\brief Deleted default constructor because this class manages an external resource.
    algorithm_executor_blocking() = delete;
    //!\brief This class provides unique ownership over the managed resource and is therefor not copyable.
    algorithm_executor_blocking(algorithm_executor_blocking const &) = delete;

    /*!\brief Move constructs the resource of the other executor.
     * \param[in] other The other executor (prvalue) to move from.
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
    algorithm_executor_blocking(algorithm_executor_blocking && other) noexcept :
        algorithm_executor_blocking{std::move(other), other.resource_position()}
    {}

    //!\brief This class provides unique ownership over the managed resource and is therefor not copyable.
    algorithm_executor_blocking & operator=(algorithm_executor_blocking const &) = delete;

    //!\brief Move assigns from the resource of another executor.
    //!\copydetails seqan3::detail::algorithm_executor_blocking::algorithm_executor_blocking(algorithm_executor_blocking && other)
    algorithm_executor_blocking & operator=(algorithm_executor_blocking && other)
    {
        auto old_resource_position = other.resource_position();

        resource = std::move(other.resource);
        move_initialise(std::move(other), old_resource_position);
        return *this;
    }

    //!\brief Defaulted.
    ~algorithm_executor_blocking() = default;

    /*!\brief Constructs this executor with the given resource range.
     *
     * \param[in] resource The underlying resource.
     * \param[in] algorithm The algorithm to invoke on the elements of the underlying resource.
     * \param[in] result A dummy result object to deduce the type of the underlying buffer value.
     * \param[in] exec_handler The execution handler to use [optional].
     *
     * \details
     *
     * If the execution handler is parallel, it allocates a buffer of the size of the given resource range.
     * Otherwise the buffer size is 1.
     * Also note that the third argument is used for deducing the algorithm result type and is otherwise
     * not used in the context of the class' construction.
     */
    algorithm_executor_blocking(resource_t resource,
                                algorithm_t algorithm,
                                algorithm_result_t const SEQAN3_DOXYGEN_ONLY(result) = algorithm_result_t{},
                                execution_handler_t && exec_handler = execution_handler_t{}) :
        exec_handler{std::move(exec_handler)},
        resource{std::views::all(resource)},
        resource_it{std::ranges::begin(this->resource)},
        algorithm{std::move(algorithm)}
    {
        if constexpr (std::same_as<execution_handler_t, execution_handler_parallel>)
            buffer_size = static_cast<size_t>(std::ranges::distance(this->resource));

        buffer.resize(buffer_size);
        buffer_it = buffer.end();
        buffer_end_it = buffer_it;
    }
    //!}

    /*!\brief Returns the next available algorithm result.
     * \returns A std::optional that either contains the next algorithm result or is empty, i.e. the
     *          underlying resource has been completely consumed.
     *
     * \details
     *
     * If there is no algorithm result available anymore the buffer will be refilled until either there is a new
     * result available or the end of the underlying resource was reached.
     * This operation is blocking such that the next result is only available after all algorithm invocations
     * triggered during seqan3::detail::algorithm_executor_blocking::fill_buffer have finished.
     *
     * ### Exception
     *
     * Throws std::bad_function_call if the algorithm was not set.
     */
    std::optional<algorithm_result_t> next_result()
    {
        fill_status status;
        // Each invocation of the algorithm might produce zero results (e.g. a search might not find a query)
        // this repeats the algorithm until it produces the first result or the input resource was consumed.
        do { status = fill_buffer(); } while (status == fill_status::empty_buffer);

        if (status == fill_status::end_of_resource)
            return {std::nullopt};

        assert(status == fill_status::non_empty_buffer);
        assert(bucket_it != buffer_it->end());

        std::optional<algorithm_result_t> result = std::ranges::iter_move(bucket_it);
        go_to_next_result(); // Go to next buffered result.
        return result;
    }

    //!\brief Checks whether the end of the input resource was reached.
    bool is_eof() noexcept
    {
        return resource_it == std::ranges::end(resource);
    }

private:
    /*!\brief This constructor is needed to ensure initialisation order for the move construction.
     * \details
     * We need to access the processed seqan3::detail::algorithm_executor_blocking::resource_position BEFORE moving the
     * resource range. As that could invalidate iterators.
     */
    algorithm_executor_blocking(algorithm_executor_blocking && other,
                                resource_difference_type old_resource_position) noexcept :
        resource{std::move(other.resource)}
    {
        move_initialise(std::move(other), old_resource_position);
    }

    /*!\brief How many times was the resource_it incremented?
     * \details
     * \note This function isn't const-qualified, as seqan3::detail::algorithm_executor_blocking::resource isn't
     * required to be seqan3::const_iterable_range by this class (i.e. having a const-qualified begin/end member
     * function). The same reason why seqan3::detail::algorithm_executor_blocking::is_eof() isn't const-qualified.
     */
    resource_difference_type resource_position()
    {
        // Get the old resource position.
        auto position = std::ranges::distance(std::ranges::begin(resource), resource_it);
        assert (position >= 0);
        return position;
    }

    //!\brief Fills the buffer by storing the results of an algorithm invocation into a pre-assigned bucket.
    fill_status fill_buffer()
    {
        if (!is_buffer_empty())  // Not everything consumed yet.
            return fill_status::non_empty_buffer;

        if (is_eof())  // Case: reached end of resource.
            return fill_status::end_of_resource;

        // Reset the buckets and the buffer iterator.
        reset_buffer();

        // Execute the algorithm (possibly asynchronous) and fill the buckets in this pre-assigned order.
        for (buffer_end_it = buffer_it; buffer_end_it != buffer.end() && !is_eof(); ++buffer_end_it, ++resource_it)
        {
            exec_handler.execute(algorithm, *resource_it, [target_buffer_it = buffer_end_it] (auto && algorithm_result)
            {
                target_buffer_it->push_back(std::move(algorithm_result));
            });
        }

        exec_handler.wait();

        // Move the results iterator to the next available result. (This skips empty results of the algorithm)
        find_next_non_empty_bucket();

        if (is_buffer_empty())
            return fill_status::empty_buffer;

        return fill_status::non_empty_buffer;
    }

    /*!\brief Whether the internal buffer is empty.
     * \returns `true` if all elements of the internal buffer have been consumed, otherwise `false`.
     */
    bool is_buffer_empty() const
    {
        return buffer_it == buffer_end_it;
    }

    /*!\brief Resets the buckets.
     *
     * \details
     *
     * Clears all buckets and sets the buffer iterator to the first bucket. The buckets are not shrunk such that the
     * allocated memory for each bucket can be reused between invocations of
     * seqan3::detail::algorithm_executor_blocking::fill_buffer.
     */
    void reset_buffer()
    {
        // Clear all buckets
        for (auto & bucket : buffer)
            bucket.clear();

        // Reset the iterator over the buckets.
        buffer_it = buffer.begin();
    }

    /*!\brief Finds the first non-empty bucket starting from the current position of the buffer iterator.
     *
     * \details
     *
     * Finds the first non-empty bucket and sets the bucket iterator to the first element of this bucket.
     * If all buckets are empty, then the buffer iterator is set to the end of the buffer and the bucket
     * iterator is not modified.
     */
    void find_next_non_empty_bucket()
    {
        assert(buffer_it <= buffer_end_it);
        // find first buffered bucket that contains at least one element
        buffer_it = std::find_if(buffer_it, buffer_end_it, [] (auto const & buffer)
        {
            return !buffer.empty();
        });

        if (buffer_it != buffer_end_it)
            bucket_it = buffer_it->begin();
    }

    /*!\brief Moves the bucket iterator to the next available result.
     *
     * \details
     *
     * If the current bucket is consumed, then the buffer iterator is incremented and the next non-empty bucket is found
     * by calling seqan3::detail::algorithm_executor_blocking::find_next_non_empty_bucket.
     */
    void go_to_next_result()
    {
        if (++bucket_it == buffer_it->end())
        {
            ++buffer_it;
            find_next_non_empty_bucket();
        }
    }

    //!\brief Helper function to move initialise `this` from `other`.
    void move_initialise(algorithm_executor_blocking && other, resource_difference_type old_resource_position) noexcept
    {
        algorithm = std::move(other.algorithm);
        buffer_size = std::move(other.buffer_size);
        exec_handler = std::move(other.exec_handler);
        // Move the resource and set the iterator state accordingly.
        resource_it = std::ranges::next(std::ranges::begin(resource), old_resource_position);

        // Get the old buffer and bucket iterator positions.
        auto buffer_it_position = other.buffer_it - other.buffer.begin();
        auto buffer_end_it_position = other.buffer_end_it - other.buffer.begin();

        std::ptrdiff_t bucket_it_position = 0;
        if (buffer_it_position != buffer_end_it_position)
            bucket_it_position = other.bucket_it - other.buffer_it->begin();

        // Move the buffer and set the buffer and bucket iterator accordingly.
        buffer = std::move(other.buffer);
        buffer_it = buffer.begin() + buffer_it_position;
        buffer_end_it = buffer.begin() + buffer_end_it_position;

        if (buffer_it_position != buffer_end_it_position)
            bucket_it = buffer_it->begin() + bucket_it_position;
    }

    //!\brief The execution policy.
    execution_handler_t exec_handler{};

    //!\brief The underlying resource.
    resource_type resource;
    //!\brief The iterator over the resource that stores the current state of the executor.
    resource_iterator_type resource_it{};
    //!\brief The algorithm to invoke.
    algorithm_t algorithm{};

    //!\brief The buffer storing the algorithm results in buckets.
    buffer_type buffer{};
    //!\brief The iterator pointing to the current bucket in the buffer.
    buffer_iterator_type buffer_it{};
    //!\brief The iterator pointing behind the last bucket (must not be the end of the buffer).
    buffer_iterator_type buffer_end_it{};
    //!\brief The bucket iterator pointing to the next result within the current bucket.
    bucket_iterator_type bucket_it{};
    //!\brief The end get pointer in the buffer.
    size_t buffer_size{1};
};

/*!\name Type deduction guides
 * \relates seqan3::detail::algorithm_executor_blocking
 * \{
 */

//!\brief Deduce the type from the provided arguments and set the sequential execution handler.
template <typename resource_rng_t, std::semiregular algorithm_t, std::semiregular algorithm_result_t>
algorithm_executor_blocking(resource_rng_t &&, algorithm_t, algorithm_result_t const &) ->
    algorithm_executor_blocking<resource_rng_t, algorithm_t, algorithm_result_t, execution_handler_sequential>;
//!\}
} // namespace seqan3::detail

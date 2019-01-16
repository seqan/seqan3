// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

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
#include <seqan3/range/shortcuts.hpp>
#include <seqan3/std/ranges>

namespace seqan3::detail
{

/*!\brief A two way executor for pairwise alignments.
 * \ingroup execution
 * \tparam align_instance_rng_t  A bounded range over the alignment instances to compute.
 * \tparam alignment_seclector_t Selects the alignment to execute.
 * \tparam execution_handler_t   The execution handler managing the execution of the alignment instances.
 */
template <std::ranges::ViewableRange align_instance_rng_t,
          typename alignment_selector_t,
          typename execution_handler_t = execution_handler_sequential>
class alignment_executor_two_way
{
    //!\brief Shortcut declaration for this class.
    using this_t = alignment_executor_two_way<align_instance_rng_t, alignment_selector_t, execution_handler_t>;

    //!\brief Grants alignment_range access to the protected/private member of this class.
    template <typename executor_t>
    friend class seqan3::alignment_range;

    /*!\name Resource
     * \{
     */
    //!\brief The resource type
    using resource_type       = std::remove_reference_t<align_instance_rng_t>;
    //!\brief The iterator over the resource.
    using resource_iterator   = std::ranges::iterator_t<resource_type>;
    //!\brief The sentinel over the resource.
    using resource_sentinel   = std::ranges::sentinel_t<resource_type>;
    //!\brief The value type of the resource.
    using resource_value_type = std::remove_reference_t<decltype(*std::declval<resource_iterator>())>;
    //!\}

    /*!\name Buffer
     * \{
     */
    //!\brief The result of invoking the alignment instance.
    using buffer_value_type = typename std::remove_reference_t<alignment_selector_t>::result_type;
    //!\brief The internal buffer.
    using buffer_type       = std::vector<buffer_value_type>;
    //!\brief The pointer type of the buffer.
    using buffer_pointer    = std::ranges::iterator_t<buffer_type>;
    //!\}
public:
    //!\brief The result type of invoking the alignment instance.
    using value_type      = buffer_value_type;
    //!\brief A reference to the alignment result.
    using reference       = std::add_lvalue_reference_t<value_type>;
    //!\brief The difference type for the buffer.
    using difference_type = typename buffer_type::difference_type;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    alignment_executor_two_way()                                               = default;
    alignment_executor_two_way(alignment_executor_two_way const &)             = default;
    alignment_executor_two_way(alignment_executor_two_way &&)                  = default;
    alignment_executor_two_way & operator=(alignment_executor_two_way const &) = default;
    alignment_executor_two_way & operator=(alignment_executor_two_way &&)      = default;
    ~alignment_executor_two_way()                                              = default;

    //!\brief Constructs this executor with the passed range of alignment instances.
    alignment_executor_two_way(align_instance_rng_t _resource,
                               alignment_selector_t _selector) :
        resource{std::forward<align_instance_rng_t>(_resource)},
        selector{std::forward<alignment_selector_t>(_selector)},
        resource_iter{begin(resource)},
        resource_end{end(resource)}
    {
        buffer.resize(1);
        setg(end(buffer), end(buffer));
    }
    //!}

    //!\brief Returns a seqan3::detail::alignment_range over the invoked alignment instances.
    alignment_range<this_t> range()
    {
        return alignment_range<this_t>{*this};
    }

protected:

    /*!\name Get area
     * \{
     */
    //!\brief Returns the remaining number of elements in the buffer, that are not read yet.
    constexpr size_t in_avail() const noexcept
    {
        return egptr - gptr;
    }

    //!\brief Returns the current alignment result in the buffer and advances the buffer get pointer.
    std::optional<std::reference_wrapper<value_type>> bump()
    {
        if (gptr == buffer_pointer{} || gptr >= egptr)
        {
            if (underflow() == eof)
            {
                return {};
            }
        }
        return {std::ref(*gptr++)};
    }

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
        setg(begin(buffer), begin(buffer) +
            std::min(static_cast<size_t>(resource_end - resource_iter), buffer.size()));

        // Apply the alignment execution.
        // TODO: Adapt for async behavior for parallel execution handler.
        std::transform(resource_iter, resource_iter + in_avail(), gptr,
            [this](auto && align_instance){
                // value_type tmp;
                auto f = selector.select(std::forward<decltype(align_instance)>(align_instance));
                value_type tmp;
                exec_handler.execute(std::move(f), tmp, [&tmp](auto && res){ tmp = std::move(res); });
                return tmp;
            });
        // Advance the resource.
        std::advance(resource_iter, in_avail());
        return in_avail();
    }
    //!\}

    /*!\name Miscellaneous
     * \{
     */
    //!\brief Checks whether the end of the input resource was reached.
    bool is_eof() const noexcept
    {
        return resource_iter == resource_end;
    }
    //!\}

private:

    //!\brief Indicates the end-of-stream.
    static constexpr size_t eof{-1};

    //!\brief The execution policy.
    execution_handler_t exec_handler{};

    //!\brief The underlying resource containing the alignment instances.
    align_instance_rng_t resource;  // view or lvalue ref
    //!\brief Selects the correct alignment to execute.
    alignment_selector_t selector;

    //!\brief Points to the current element in the resource.
    resource_iterator resource_iter{};
    //!\brief End of the resource.
    resource_sentinel resource_end{};

    //!\brief The buffer storing the alignment results.
    buffer_type buffer;
    //!\brief The get pointer in the buffer.
    buffer_pointer gptr;
    //!\brief The end get pointer in the buffer.
    buffer_pointer egptr;
};

/*!\name Type deduction guides
 * \relates seqan3::detail::alignment_executor_two_way
 * \{
 */
template <std::ranges::ViewableRange resource_rng_t, typename selector_t>
alignment_executor_two_way(resource_rng_t && , selector_t &&) ->
    alignment_executor_two_way<resource_rng_t, selector_t, execution_handler_sequential>;

//!\}
} // namespace seqan3::detail

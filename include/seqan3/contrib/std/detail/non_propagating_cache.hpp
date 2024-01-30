// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 * \brief Provides [seqan::stl::detail::non_propagating_cache](https://eel.is/c++draft/range.nonprop.cache).
 */

// File might be included from multiple libraries.
#ifndef SEQAN_STD_DETAIL_NON_PROPAGATING_CACHE
#define SEQAN_STD_DETAIL_NON_PROPAGATING_CACHE

#include <optional>

namespace seqan::stl::detail
{

/*!\brief A helper that enables an input view to temporarily cache values as it is iterated over.
 * \sa https://eel.is/c++draft/range.nonprop.cache
 * \details
 *
 *  Behaves like std::optional, but
 *  - Constrains its type parameter `T` with `std::is_object_v<T>`.
 *  - Copy constructor is a no-op.
 *  - Move constructor is a no-op for `this` and destroys value of `other`.
 *  - Copy assignment is a no-op if `this` == `other`, otherwise destroys value of `this`.
 *  - Move assignment destroys both objects.
 *  - Has an additional `emplace_deref` member function.
 */
template <typename T>
    requires std::is_object_v<T>
class non_propagating_cache : public std::optional<T>
{
public:
    //!\brief Use std::optional constructors.
    using std::optional<T>::optional;
    //!\brief Use std::optional assignment operators.
    using std::optional<T>::operator=;

    //!\brief Copy construction is a no-op.
    constexpr non_propagating_cache(non_propagating_cache const &) noexcept
    {}

    //!\brief Move constructor is a no-op for `this` and destroys value of `other`.
    constexpr non_propagating_cache(non_propagating_cache && other) noexcept
    {
        other.reset();
    }

    //!\brief Copy assignment is a no-op if `this` == `other`, otherwise destroys value of `this`.
    constexpr non_propagating_cache & operator=(non_propagating_cache const & other) noexcept
    {
        if (std::addressof(other) != this)
            this->reset();
        return *this;
    }

    //!\brief Move assignment destroys values of both `this` and `other`.
    constexpr non_propagating_cache & operator=(non_propagating_cache && other) noexcept
    {
        this->reset();
        other.reset();
        return *this;
    }

    //!\brief Destroys current value, initializes with new value obtained from dereferencing `i`, and returns value.
    template <class I>
    constexpr T & emplace_deref(I const & i)
        requires requires (I const & i) { T(*i); }
    {
        this->reset();
        this->emplace(*i);
        return this->value();
    }
};

} // namespace seqan::stl::detail

#endif // SEQAN_STD_DETAIL_NON_PROPAGATING_CACHE

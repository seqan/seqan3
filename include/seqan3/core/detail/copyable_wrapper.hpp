// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::detail::copyable_wrapper.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <functional>
#include <optional>

#include <seqan3/core/platform.hpp>

namespace seqan3::detail
{

/*!\brief Helper concept for seqan3::detail::copyable_wrapper.
 * \ingroup core
 * \see https://en.cppreference.com/w/cpp/ranges/copyable_wrapper
 * \noapi{Exposition only.}
 */
template <typename t>
concept boxable = std::copy_constructible<t> && std::is_object_v<t>;

/*!\brief Utility wrapper that behaves like std::optional but makes the type conform with the std::copyable concept.
 * \ingroup core
 * \see https://en.cppreference.com/w/cpp/ranges/copyable_wrapper
 */
template <typename t>
class copyable_wrapper : public std::optional<t>
{
public:
    using std::optional<t>::optional;                               //!< Use std::optional constructors.
    using std::optional<t>::operator=;                              //!< Use std::optional assignment operators.
    constexpr copyable_wrapper(copyable_wrapper const &) = default; //!< Defaulted.
    constexpr copyable_wrapper(copyable_wrapper &&) = default;      //!< Defaulted.
    constexpr ~copyable_wrapper() = default;                        //!< Defaulted.

    /*!\brief Use a specialised default constructor, if the wrapped type is default initialisable.
     *        If not, the default constructor of std::optional is used.
     */
    constexpr copyable_wrapper() noexcept(std::is_nothrow_default_constructible_v<t>)
        requires std::default_initializable<t>
        : copyable_wrapper{std::in_place}
    {}

    //!\brief Copy assignment for non-copyable wrapped types.
    constexpr copyable_wrapper & operator=(copyable_wrapper const & other)
        noexcept(std::is_nothrow_copy_constructible_v<t>)
        requires (!std::copyable<t>)
    {
        if (this != std::addressof(other))
        {
            if (other)
                this->emplace(*other);
            else
                this->reset();
        }
        return *this;
    }

    //!\brief Move assignment for non-movable wrapped types.
    constexpr copyable_wrapper & operator=(copyable_wrapper && other) noexcept(std::is_nothrow_move_constructible_v<t>)
        requires (!std::movable<t>)
    {
        if (this != std::addressof(other))
        {
            if (other)
                this->emplace(std::move(*other));
            else
                this->reset();
        }
        return *this;
    }

    /*!\brief Invokes the wrapped object with passed arguments.
     * \throws std::bad_optional_access if no value is contained.
     * \details
     *
     * This is a SeqAn-specific extension that allows easy invocation of the wrapped object.
     * `t` needs to be callable with the passed arguments.
     *
     * Constness of this function depends on the constness of the call-operator of the wrapped object.
     *
     * \include test/snippet/core/detail/copyable_wrapper.cpp
     */
    template <typename... args_t>
        requires std::invocable<t, args_t...>
    constexpr decltype(auto) operator()(args_t... args) noexcept(std::is_nothrow_invocable_v<t, args_t...>)
    {
        return std::invoke(this->value(), std::forward<args_t>(args)...);
    }

    //!\overload
    template <typename... args_t>
        requires std::invocable<t, args_t...>
    constexpr decltype(auto) operator()(args_t... args) const noexcept(std::is_nothrow_invocable_v<t, args_t...>)
    {
        return std::invoke(this->value(), std::forward<args_t>(args)...);
    }
};

/*!\brief Utility wrapper that behaves like std::optional but makes the type conform with the std::copyable concept.
 * \ingroup core
 * \see https://en.cppreference.com/w/cpp/ranges/copyable_wrapper
 * \details
 *
 * If `t` is `std::copyable`, the STL allows to store `t` directly, i.e. no `std::optional` is needed.
 */
template <boxable t>
    requires std::copyable<t>
          || (std::is_nothrow_move_constructible_v<t> && std::is_nothrow_copy_constructible_v<t>)
#if SEQAN3_WORKAROUND_DEFAULT_CONSTRUCTIBLE_VIEW // If views must be default constructible, t must also be
                 && std::default_initializable<t>
#endif
             class copyable_wrapper<t>
{
private:
    t value{}; //!< An object of the wrapped type.

public:
    //!\brief copyable_wrapper is default constructible, iff the wrapped type is default initialisable.
    constexpr copyable_wrapper()
        requires std::default_initializable<t>
    = default;

    constexpr copyable_wrapper(copyable_wrapper const &) = default; //!< Defaulted.
    constexpr copyable_wrapper(copyable_wrapper &&) = default;      //!< Defaulted.
    constexpr ~copyable_wrapper() = default;                        //!< Defaulted.

    //!\brief Copy assignment for copyable types is the default copy assignment.
    constexpr copyable_wrapper & operator=(copyable_wrapper const &)
        requires std::copyable<t>
    = default;

    //!\brief Copy assignment for non-copyable types uses the Destroy-then-copy paradigm.
    constexpr copyable_wrapper & operator=(copyable_wrapper const & other) noexcept
    {
        // Destroy-then-copy
        if (this != std::addressof(other))
        {
            value.~t();                                       // Call destructor of value
            std::construct_at(std::addressof(value), *other); // Copy construct
        }
        return *this;
    }

    //!\brief Move assignment for copyable types is the default move assignment.
    constexpr copyable_wrapper & operator=(copyable_wrapper &&)
        requires std::copyable<t>
    = default;

    //!\brief Move assignment for non-copyable types uses the Destroy-then-copy paradigm.
    constexpr copyable_wrapper & operator=(copyable_wrapper && other) noexcept
    {
        // Destroy-then-copy
        if (this != std::addressof(other))
        {
            value.~t();                                                  // Call destructor of value
            std::construct_at(std::addressof(value), std::move(*other)); // Move construct
        }
        return *this;
    }

    //!\brief Copy construct from value.
    constexpr explicit copyable_wrapper(t const & other) noexcept(std::is_nothrow_copy_constructible_v<t>) :
        value{other}
    {}

    //!\brief Move construct from value.
    constexpr explicit copyable_wrapper(t && other) noexcept(std::is_nothrow_move_constructible_v<t>) :
        value{std::move(other)}
    {}

    //!\brief Construct from argument pack. Part of the std::optional API.
    template <typename... args_t>
        requires std::constructible_from<t, args_t...>
    constexpr explicit copyable_wrapper(std::in_place_t, args_t... args)
        noexcept(std::is_nothrow_constructible_v<t, args_t...>) :
        value{std::forward<args_t>(args)...}
    {}

    //!\brief Return whether the wrapper contains a value. Part of the std::optional API.
    constexpr bool has_value() const noexcept
    {
        return true; // t is copyable, hence we always store a value.
    }

    //!\brief Returns a reference to the wrapped object. Part of the std::optional API.
    constexpr t & operator*() noexcept
    {
        return value;
    }

    //!\brief Returns a reference to the wrapped object. Part of the std::optional API.
    constexpr t const & operator*() const noexcept
    {
        return value;
    }

    //!\brief Returns a pointer to the wrapped object. Part of the std::optional API.
    constexpr t * operator->() noexcept
    {
        return std::addressof(value);
    }

    //!\brief Returns a pointer to the wrapped object. Part of the std::optional API.
    constexpr t const * operator->() const noexcept
    {
        return std::addressof(value);
    }

    /*!\brief Invokes the wrapped object with passed arguments.
     * \details
     *
     * This is a SeqAn-specific extension that allows easy invocation of the wrapped object.
     * `t` needs to be callable with the passed arguments.
     *
     * Constness of this function depends on the constness of the call-operator of the wrapped object.
     *
     * \include test/snippet/core/detail/copyable_wrapper.cpp
     */
    template <typename... args_t>
        requires std::invocable<t, args_t...>
    constexpr decltype(auto) operator()(args_t... args) noexcept(std::is_nothrow_invocable_v<t, args_t...>)
    {
        return std::invoke(value, std::forward<args_t>(args)...);
    }

    //!\overload
    template <typename... args_t>
        requires std::invocable<t, args_t...>
    constexpr decltype(auto) operator()(args_t... args) const noexcept(std::is_nothrow_invocable_v<t, args_t...>)
    {
        return std::invoke(value, std::forward<args_t>(args)...);
    }
};

/*!\brief Type deduction guide that strips references.
 * \relates seqan3::detail::copyable_wrapper
 */
template <typename t>
copyable_wrapper(t) -> copyable_wrapper<std::remove_reference_t<t>>;

/*!\brief Utility transformation trait to get a wrapper type that models std::copyable.
 * \ingroup core
 *
 * \see https://en.cppreference.com/w/cpp/ranges/copyable_wrapper
 */
template <typename t>
using copyable_wrapper_t = copyable_wrapper<std::remove_reference_t<t>>;

} // namespace seqan3::detail

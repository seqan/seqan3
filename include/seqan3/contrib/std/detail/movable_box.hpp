// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 * \brief Provides seqan::stl::detail::movable_box.
 */

// File might be included from multiple libraries.
#ifndef SEQAN_STD_DETAIL_MOVABLE_BOX
#define SEQAN_STD_DETAIL_MOVABLE_BOX

#include <functional>
#include <optional>

namespace seqan::stl::detail
{

template <typename t>
concept boxable = std::move_constructible<t> && std::is_object_v<t>;

template <boxable t>
class movable_box : public std::optional<t>
{
public:
    using std::optional<t>::optional;
    using std::optional<t>::operator=;
    constexpr movable_box(movable_box const &) = default;
    constexpr movable_box(movable_box &&) = default;
    constexpr ~movable_box() = default;

    constexpr movable_box() noexcept(std::is_nothrow_default_constructible_v<t>)
        requires std::default_initializable<t>
        : movable_box{std::in_place}
    {}

    constexpr movable_box & operator=(movable_box const & other) noexcept(std::is_nothrow_copy_constructible_v<t>)
        requires (!std::copyable<t> && std::copy_constructible<t>)
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

    constexpr movable_box & operator=(movable_box && other) noexcept(std::is_nothrow_move_constructible_v<t>)
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

    template <typename... args_t>
        requires std::invocable<t, args_t...>
    constexpr decltype(auto) operator()(args_t... args) noexcept(std::is_nothrow_invocable_v<t, args_t...>)
    {
        return std::invoke(this->value(), std::forward<args_t>(args)...);
    }

    template <typename... args_t>
        requires std::invocable<t, args_t...>
    constexpr decltype(auto) operator()(args_t... args) const noexcept(std::is_nothrow_invocable_v<t, args_t...>)
    {
        return std::invoke(this->value(), std::forward<args_t>(args)...);
    }
};

template <boxable t>
    requires std::copyable<t> || (std::is_nothrow_move_constructible_v<t> && std::is_nothrow_copy_constructible_v<t>)
class movable_box<t>
{
private:
    t value{};

public:
    constexpr movable_box()
        requires std::default_initializable<t>
    = default;

    constexpr movable_box(movable_box const &) = default;
    constexpr movable_box(movable_box &&) = default;
    constexpr ~movable_box() = default;

    constexpr movable_box & operator=(movable_box const &) = default;

    constexpr movable_box & operator=(movable_box &&) = default;

    constexpr explicit movable_box(t const & other) noexcept(std::is_nothrow_copy_constructible_v<t>) : value{other}
    {}

    constexpr explicit movable_box(t && other) noexcept(std::is_nothrow_move_constructible_v<t>) :
        value{std::move(other)}
    {}

    template <typename... args_t>
        requires std::constructible_from<t, args_t...>
    constexpr explicit movable_box(std::in_place_t,
                                   args_t... args) noexcept(std::is_nothrow_constructible_v<t, args_t...>) :
        value{std::forward<args_t>(args)...}
    {}

    constexpr bool has_value() const noexcept
    {
        return true; // t is copyable, hence we always store a value.
    }

    constexpr t & operator*() noexcept
    {
        return value;
    }

    constexpr t const & operator*() const noexcept
    {
        return value;
    }

    constexpr t * operator->() noexcept
    {
        return std::addressof(value);
    }

    constexpr t const * operator->() const noexcept
    {
        return std::addressof(value);
    }

    template <typename... args_t>
        requires std::invocable<t, args_t...>
    constexpr decltype(auto) operator()(args_t... args) noexcept(std::is_nothrow_invocable_v<t, args_t...>)
    {
        return std::invoke(value, std::forward<args_t>(args)...);
    }

    template <typename... args_t>
        requires std::invocable<t, args_t...>
    constexpr decltype(auto) operator()(args_t... args) const noexcept(std::is_nothrow_invocable_v<t, args_t...>)
    {
        return std::invoke(value, std::forward<args_t>(args)...);
    }
};

template <typename t>
movable_box(t) -> movable_box<std::remove_reference_t<t>>;

template <typename t>
using movable_box_t = movable_box<std::remove_reference_t<t>>;

} // namespace seqan::stl::detail

#endif // SEQAN_STD_DETAIL_MOVABLE_BOX

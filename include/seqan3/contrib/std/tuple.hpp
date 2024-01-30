// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan::stl::tuple.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#ifndef SEQAN_STD_TUPLE
#define SEQAN_STD_TUPLE

#include <tuple>

#ifdef __cpp_lib_tuple_like

namespace seqan::stl
{

using std::tuple;

} // namespace seqan::stl

#else

#    include "detail/exposition_only.hpp"
#    include "pair.hpp"

namespace seqan::stl
{

template <class... Types>
class tuple : public std::tuple<Types...>
{
private:
    //!\brief The underlying std::tuple type.
    using base_t = std::tuple<Types...>;

    //!\brief Constructs by unfolding another tuple-like object.
    template <typename tuple_like_t, std::size_t... N>
    tuple(tuple_like_t && other, std::integer_sequence<size_t, N...>) :
        base_t((std::forward<std::tuple_element_t<N, tuple_like_t>>(std::get<N>(other)))...)
    {}

public:
    /*!\name Default constructors and assignments.
     * \{
     */
    tuple() = default;                          //!< Defaulted.
    tuple(tuple const &) = default;             //!< Defaulted.
    tuple & operator=(tuple const &) = default; //!< Defaulted.
    ~tuple() = default;                         //!< Defaulted.
    //!\}

    /*!\brief Returns the tuple as the underlying std::tuple type.
     * \private
     */
    base_t & as_base() noexcept
    {
        return *this;
    }

    /*!\brief Returns the tuple as the underlying std::tuple type.
     * \private
     */
    base_t const & as_base() const noexcept
    {
        return *this;
    }

    /*!\name Construct from arguments.
     * \{
     */
    //!@{ Constructs from arguments.
    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_constructible_v<Types, UTypes &> && ...)
    constexpr tuple(UTypes &... other) : base_t(other...)
    {}

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_constructible_v<Types, UTypes const> && ...)
    constexpr tuple(UTypes const &... other) : base_t(other...)
    {}

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_constructible_v<Types, UTypes> && ...)
    constexpr tuple(UTypes &&... other) : base_t(std::forward<UTypes>(other)...)
    {}

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_constructible_v<Types, UTypes const> && ...)
    constexpr tuple(UTypes const &&... other) : base_t(std::forward<UTypes const>(other)...)
    {}
    //!@}
    //!\}

    /*!\name Construct from tuple.
     * \{
     */
    //!@{ Constructs from tuple.
    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_constructible_v<Types, UTypes &> && ...)
    constexpr tuple(tuple<UTypes...> & other) : tuple(other, std::index_sequence_for<Types...>{})
    {}

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_constructible_v<Types, UTypes const> && ...)
    constexpr tuple(tuple<UTypes...> const & other) : tuple(other, std::index_sequence_for<Types...>{})
    {}

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_constructible_v<Types, UTypes> && ...)
    constexpr tuple(tuple<UTypes...> && other) : tuple(std::move(other), std::index_sequence_for<Types...>{})
    {}

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_constructible_v<Types, UTypes const> && ...)
    constexpr tuple(tuple<UTypes...> const && other) : tuple(std::move(other), std::index_sequence_for<Types...>{})
    {}
    //!@}
    //!\}

    /*!\name Construct from pair.
     * \{
     */
    //!@{ Constructs from pair.
    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_constructible_v<Types, UTypes &> && ...)
    constexpr tuple(pair<UTypes...> & other) : tuple(other, std::index_sequence_for<Types...>{})
    {}

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_constructible_v<Types, UTypes const> && ...)
    constexpr tuple(pair<UTypes...> const & other) : tuple(other, std::index_sequence_for<Types...>{})
    {}

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_constructible_v<Types, UTypes> && ...)
    constexpr tuple(pair<UTypes...> && other) : tuple(std::move(other), std::index_sequence_for<Types...>{})
    {}

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_constructible_v<Types, UTypes const> && ...)
    constexpr tuple(pair<UTypes...> const && other) : tuple(std::move(other), std::index_sequence_for<Types...>{})
    {}
    //!@}
    //!\}

    /*!\name Construct from std::tuple.
     * \{
     */
    //!@{ Constructs from std::tuple.
    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_constructible_v<Types, UTypes &> && ...)
    constexpr tuple(std::tuple<UTypes...> & other) : tuple(other, std::index_sequence_for<Types...>{})
    {}

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_constructible_v<Types, UTypes const> && ...)
    constexpr tuple(std::tuple<UTypes...> const & other) : tuple(other, std::index_sequence_for<Types...>{})
    {}

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_constructible_v<Types, UTypes> && ...)
    constexpr tuple(std::tuple<UTypes...> && other) : tuple(std::move(other), std::index_sequence_for<Types...>{})
    {}

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_constructible_v<Types, UTypes const> && ...)
    constexpr tuple(std::tuple<UTypes...> const && other) : tuple(std::move(other), std::index_sequence_for<Types...>{})
    {}

    /*!\name Construct from std::pair.
     * \{
     */
    //!@{ Constructs from std::pair.
    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_constructible_v<Types, UTypes &> && ...)
    constexpr tuple(std::pair<UTypes...> & other) : tuple(other, std::index_sequence_for<Types...>{})
    {}

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_constructible_v<Types, UTypes const> && ...)
    constexpr tuple(std::pair<UTypes...> const & other) : tuple(other, std::index_sequence_for<Types...>{})
    {}

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_constructible_v<Types, UTypes> && ...)
    constexpr tuple(std::pair<UTypes...> && other) : tuple(std::move(other), std::index_sequence_for<Types...>{})
    {}

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_constructible_v<Types, UTypes const> && ...)
    constexpr tuple(std::pair<UTypes...> const && other) : tuple(std::move(other), std::index_sequence_for<Types...>{})
    {}
    //!@}
    //!\}

    /*!\name Assign from tuple.
     * \{
     */
    //!@{ Assigns from tuple.
    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_assignable_v<Types, UTypes &> && ...)
    constexpr tuple & operator=(tuple<UTypes...> & other)
    {
        [&]<size_t... N>(std::integer_sequence<size_t, N...>)
        {
            ((std::get<N>(*this) = std::get<N>(other)), ...);
        }
        (std::index_sequence_for<Types...>{});

        return *this;
    }

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_assignable_v<Types, UTypes const> && ...)
    constexpr tuple & operator=(tuple<UTypes...> const & other)
    {
        [&]<size_t... N>(std::integer_sequence<size_t, N...>)
        {
            ((std::get<N>(*this) = std::get<N>(other)), ...);
        }
        (std::index_sequence_for<Types...>{});

        return *this;
    }

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_assignable_v<Types, UTypes> && ...)
    constexpr tuple & operator=(tuple<UTypes...> && other)
    {
        [&]<size_t... N>(std::integer_sequence<size_t, N...>)
        {
            ((std::get<N>(*this) = std::get<N>(other)), ...);
        }
        (std::index_sequence_for<Types...>{});

        return *this;
    }

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_assignable_v<Types, UTypes const> && ...)
    constexpr tuple & operator=(tuple<UTypes...> const && other)
    {
        [&]<size_t... N>(std::integer_sequence<size_t, N...>)
        {
            ((std::get<N>(*this) = std::get<N>(other)), ...);
        }
        (std::index_sequence_for<Types...>{});

        return *this;
    }

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_assignable_v<Types const, UTypes &> && ...)
    constexpr tuple const & operator=(tuple<UTypes...> & other) const
    {
        [&]<size_t... N>(std::integer_sequence<size_t, N...>)
        {
            ((std::get<N>(*this) = std::get<N>(other)), ...);
        }
        (std::index_sequence_for<Types...>{});

        return *this;
    }

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_assignable_v<Types const, UTypes const> && ...)
    constexpr tuple const & operator=(tuple<UTypes...> const & other) const
    {
        [&]<size_t... N>(std::integer_sequence<size_t, N...>)
        {
            ((std::get<N>(*this) = std::get<N>(other)), ...);
        }
        (std::index_sequence_for<Types...>{});

        return *this;
    }

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_assignable_v<Types const, UTypes> && ...)
    constexpr tuple const & operator=(tuple<UTypes...> && other) const
    {
        [&]<size_t... N>(std::integer_sequence<size_t, N...>)
        {
            ((std::get<N>(*this) = std::get<N>(other)), ...);
        }
        (std::index_sequence_for<Types...>{});

        return *this;
    }

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_assignable_v<Types const, UTypes const> && ...)
    constexpr tuple const & operator=(tuple<UTypes...> const && other) const
    {
        [&]<size_t... N>(std::integer_sequence<size_t, N...>)
        {
            ((std::get<N>(*this) = std::get<N>(other)), ...);
        }
        (std::index_sequence_for<Types...>{});

        return *this;
    }
    //!@}
    //!\}

    /*!\name Assign from pair.
     * \{
     */
    //!@{ Assigns from pair.
    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_assignable_v<Types, UTypes &> && ...)
    constexpr tuple & operator=(pair<UTypes...> & other)
    {
        [&]<size_t... N>(std::integer_sequence<size_t, N...>)
        {
            ((std::get<N>(*this) = std::get<N>(other)), ...);
        }
        (std::index_sequence_for<Types...>{});

        return *this;
    }

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_assignable_v<Types, UTypes const> && ...)
    constexpr tuple & operator=(pair<UTypes...> const & other)
    {
        [&]<size_t... N>(std::integer_sequence<size_t, N...>)
        {
            ((std::get<N>(*this) = std::get<N>(other)), ...);
        }
        (std::index_sequence_for<Types...>{});

        return *this;
    }

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_assignable_v<Types, UTypes> && ...)
    constexpr tuple & operator=(pair<UTypes...> && other)
    {
        [&]<size_t... N>(std::integer_sequence<size_t, N...>)
        {
            ((std::get<N>(*this) = std::get<N>(other)), ...);
        }
        (std::index_sequence_for<Types...>{});

        return *this;
    }

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_assignable_v<Types, UTypes const> && ...)
    constexpr tuple & operator=(pair<UTypes...> const && other)
    {
        [&]<size_t... N>(std::integer_sequence<size_t, N...>)
        {
            ((std::get<N>(*this) = std::get<N>(other)), ...);
        }
        (std::index_sequence_for<Types...>{});

        return *this;
    }

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_assignable_v<Types const, UTypes &> && ...)
    constexpr tuple const & operator=(pair<UTypes...> & other) const
    {
        [&]<size_t... N>(std::integer_sequence<size_t, N...>)
        {
            ((std::get<N>(*this) = std::get<N>(other)), ...);
        }
        (std::index_sequence_for<Types...>{});

        return *this;
    }

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_assignable_v<Types const, UTypes const> && ...)
    constexpr tuple const & operator=(pair<UTypes...> const & other) const
    {
        [&]<size_t... N>(std::integer_sequence<size_t, N...>)
        {
            ((std::get<N>(*this) = std::get<N>(other)), ...);
        }
        (std::index_sequence_for<Types...>{});

        return *this;
    }

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_assignable_v<Types const, UTypes> && ...)
    constexpr tuple const & operator=(pair<UTypes...> && other) const
    {
        [&]<size_t... N>(std::integer_sequence<size_t, N...>)
        {
            ((std::get<N>(*this) = std::get<N>(other)), ...);
        }
        (std::index_sequence_for<Types...>{});

        return *this;
    }

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_assignable_v<Types const, UTypes const> && ...)
    constexpr tuple const & operator=(pair<UTypes...> const && other) const
    {
        [&]<size_t... N>(std::integer_sequence<size_t, N...>)
        {
            ((std::get<N>(*this) = std::get<N>(other)), ...);
        }
        (std::index_sequence_for<Types...>{});

        return *this;
    }
    //!@}
    //!\}

    /*!\name Assign from std::tuple.
     * \{
     */
    //!@{ Assigns from std::tuple.
    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_assignable_v<Types, UTypes &> && ...)
    constexpr tuple & operator=(std::tuple<UTypes...> & other)
    {
        [&]<size_t... N>(std::integer_sequence<size_t, N...>)
        {
            ((std::get<N>(*this) = std::get<N>(other)), ...);
        }
        (std::index_sequence_for<Types...>{});

        return *this;
    }

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_assignable_v<Types, UTypes const> && ...)
    constexpr tuple & operator=(std::tuple<UTypes...> const & other)
    {
        [&]<size_t... N>(std::integer_sequence<size_t, N...>)
        {
            ((std::get<N>(*this) = std::get<N>(other)), ...);
        }
        (std::index_sequence_for<Types...>{});

        return *this;
    }

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_assignable_v<Types, UTypes> && ...)
    constexpr tuple & operator=(std::tuple<UTypes...> && other)
    {
        [&]<size_t... N>(std::integer_sequence<size_t, N...>)
        {
            ((std::get<N>(*this) = std::get<N>(other)), ...);
        }
        (std::index_sequence_for<Types...>{});

        return *this;
    }

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_assignable_v<Types, UTypes const> && ...)
    constexpr tuple & operator=(std::tuple<UTypes...> const && other)
    {
        [&]<size_t... N>(std::integer_sequence<size_t, N...>)
        {
            ((std::get<N>(*this) = std::get<N>(other)), ...);
        }
        (std::index_sequence_for<Types...>{});

        return *this;
    }

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_assignable_v<Types const, UTypes &> && ...)
    constexpr tuple const & operator=(std::tuple<UTypes...> & other) const
    {
        [&]<size_t... N>(std::integer_sequence<size_t, N...>)
        {
            ((std::get<N>(*this) = std::get<N>(other)), ...);
        }
        (std::index_sequence_for<Types...>{});

        return *this;
    }

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_assignable_v<Types const, UTypes const> && ...)
    constexpr tuple const & operator=(std::tuple<UTypes...> const & other) const
    {
        [&]<size_t... N>(std::integer_sequence<size_t, N...>)
        {
            ((std::get<N>(*this) = std::get<N>(other)), ...);
        }
        (std::index_sequence_for<Types...>{});

        return *this;
    }

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_assignable_v<Types const, UTypes> && ...)
    constexpr tuple const & operator=(std::tuple<UTypes...> && other) const
    {
        [&]<size_t... N>(std::integer_sequence<size_t, N...>)
        {
            ((std::get<N>(*this) = std::get<N>(other)), ...);
        }
        (std::index_sequence_for<Types...>{});

        return *this;
    }

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_assignable_v<Types const, UTypes const> && ...)
    constexpr tuple const & operator=(std::tuple<UTypes...> const && other) const
    {
        [&]<size_t... N>(std::integer_sequence<size_t, N...>)
        {
            ((std::get<N>(*this) = std::get<N>(other)), ...);
        }
        (std::index_sequence_for<Types...>{});

        return *this;
    }
    //!@}
    //!\}

    /*!\name Assign from std::pair.
     * \{
     */
    //!@{ Assigns from std::pair.
    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_assignable_v<Types, UTypes &> && ...)
    constexpr tuple & operator=(std::pair<UTypes...> & other)
    {
        [&]<size_t... N>(std::integer_sequence<size_t, N...>)
        {
            ((std::get<N>(*this) = std::get<N>(other)), ...);
        }
        (std::index_sequence_for<Types...>{});

        return *this;
    }

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_assignable_v<Types, UTypes const> && ...)
    constexpr tuple & operator=(std::pair<UTypes...> const & other)
    {
        [&]<size_t... N>(std::integer_sequence<size_t, N...>)
        {
            ((std::get<N>(*this) = std::get<N>(other)), ...);
        }
        (std::index_sequence_for<Types...>{});

        return *this;
    }

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_assignable_v<Types, UTypes> && ...)
    constexpr tuple & operator=(std::pair<UTypes...> && other)
    {
        [&]<size_t... N>(std::integer_sequence<size_t, N...>)
        {
            ((std::get<N>(*this) = std::get<N>(other)), ...);
        }
        (std::index_sequence_for<Types...>{});

        return *this;
    }

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_assignable_v<Types, UTypes const> && ...)
    constexpr tuple & operator=(std::pair<UTypes...> const && other)
    {
        [&]<size_t... N>(std::integer_sequence<size_t, N...>)
        {
            ((std::get<N>(*this) = std::get<N>(other)), ...);
        }
        (std::index_sequence_for<Types...>{});

        return *this;
    }

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_assignable_v<Types const, UTypes &> && ...)
    constexpr tuple const & operator=(std::pair<UTypes...> & other) const
    {
        [&]<size_t... N>(std::integer_sequence<size_t, N...>)
        {
            ((std::get<N>(*this) = std::get<N>(other)), ...);
        }
        (std::index_sequence_for<Types...>{});

        return *this;
    }

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_assignable_v<Types const, UTypes const> && ...)
    constexpr tuple const & operator=(std::pair<UTypes...> const & other) const
    {
        [&]<size_t... N>(std::integer_sequence<size_t, N...>)
        {
            ((std::get<N>(*this) = std::get<N>(other)), ...);
        }
        (std::index_sequence_for<Types...>{});

        return *this;
    }

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_assignable_v<Types const, UTypes> && ...)
    constexpr tuple const & operator=(std::pair<UTypes...> && other) const
    {
        [&]<size_t... N>(std::integer_sequence<size_t, N...>)
        {
            ((std::get<N>(*this) = std::get<N>(other)), ...);
        }
        (std::index_sequence_for<Types...>{});

        return *this;
    }

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_assignable_v<Types const, UTypes const> && ...)
    constexpr tuple const & operator=(std::pair<UTypes...> const && other) const
    {
        [&]<size_t... N>(std::integer_sequence<size_t, N...>)
        {
            ((std::get<N>(*this) = std::get<N>(other)), ...);
        }
        (std::index_sequence_for<Types...>{});

        return *this;
    }
    //!@}
    //!\}

    /*!\name Conversion to std::tuple.
     * \{
     */
    //!@{ Converts to std::tuple
    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_constructible_v<UTypes, Types &> && ...)
    operator std::tuple<UTypes...>() &
    {
        return std::make_from_tuple<std::tuple<UTypes...>>(*this);
    }

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_constructible_v<UTypes, Types const> && ...)
    operator std::tuple<UTypes...>() const &
    {
        return std::make_from_tuple<std::tuple<UTypes...>>(*this);
    }

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_constructible_v<UTypes, Types> && ...)
    operator std::tuple<UTypes...>() &&
    {
        return std::make_from_tuple<std::tuple<UTypes...>>(std::move(*this));
    }

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_constructible_v<UTypes, Types const> && ...)
    operator std::tuple<UTypes...>() const &&
    {
        return std::make_from_tuple<std::tuple<UTypes...>>(std::move(*this));
    }
    //!@}
    //!\}

    /*!\name Comparison operators (tuple)
     * \{
     */
    /*!\brief Checks whether `lhs` and `rhs` are equal.
     * \tparam UTypes The types of the elements of `rhs`. Automatically deduced.
     * \param lhs A tuple.
     * \param rhs A tuple with possibly different element types.
     * \returns A bool indicating the result of the comparison.
     */
    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes))
              && requires { requires (std::equality_comparable_with<Types, UTypes> && ...); } // Defer instantiation
    constexpr friend bool operator==(tuple const & lhs, tuple<UTypes...> const & rhs)
    {
        static_assert((std::equality_comparable_with<Types, UTypes> && ...));
        return lhs.as_base() == rhs.as_base();
    }

    /*!\brief Checks whether `lhs` and `rhs` are unequal.
     * \tparam UTypes The types of the elements of `rhs`. Automatically deduced.
     * \param lhs A tuple.
     * \param rhs A tuple with possibly different element types.
     * \returns A bool indicating the result of the comparison.
     */
    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes))
              && requires { requires (std::equality_comparable_with<Types, UTypes> && ...); } // Defer instantiation
    constexpr friend bool operator!=(tuple const & lhs, tuple<UTypes...> const & rhs)
    {
        return lhs.as_base() != rhs.as_base();
    }

    /*!\brief Checks whether `lhs` is less than `rhs`.
     * \tparam UTypes The types of the elements of `rhs`. Automatically deduced.
     * \param lhs A tuple.
     * \param rhs A tuple with possibly different element types.
     * \returns A bool indicating the result of the comparison.
     */
    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes))
              && requires { requires (std::totally_ordered_with<Types, UTypes> && ...); } // Defer instantiation
    constexpr friend bool operator<(tuple const & lhs, tuple<UTypes...> const & rhs)
    {
        return lhs.as_base() < rhs.as_base();
    }

    /*!\brief Checks whether `lhs` is less than or equal to `rhs`.
     * \tparam UTypes The types of the elements of `rhs`. Automatically deduced.
     * \param lhs A tuple.
     * \param rhs A tuple with possibly different element types.
     * \returns A bool indicating the result of the comparison.
     */
    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes))
              && requires { requires (std::totally_ordered_with<Types, UTypes> && ...); } // Defer instantiation
    constexpr friend bool operator<=(tuple const & lhs, tuple<UTypes...> const & rhs)
    {
        return lhs.as_base() <= rhs.as_base();
    }

    /*!\brief Checks whether `lhs` is greater than `rhs`.
     * \tparam UTypes The types of the elements of `rhs`. Automatically deduced.
     * \param lhs A tuple.
     * \param rhs A tuple with possibly different element types.
     * \returns A bool indicating the result of the comparison.
     */
    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes))
              && requires { requires (std::totally_ordered_with<Types, UTypes> && ...); } // Defer instantiation
    constexpr friend bool operator>(tuple const & lhs, tuple<UTypes...> const & rhs)
    {
        return lhs.as_base() > rhs.as_base();
    }

    /*!\brief Checks whether `lhs` is greater than or equal to `rhs`.
     * \tparam UTypes The types of the elements of `rhs`. Automatically deduced.
     * \param lhs A tuple.
     * \param rhs A tuple with possibly different element types.
     * \returns A bool indicating the result of the comparison.
     */
    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes))
              && requires { requires (std::totally_ordered_with<Types, UTypes> && ...); } // Defer instantiation
    constexpr friend bool operator>=(tuple const & lhs, tuple<UTypes...> const & rhs)
    {
        return lhs.as_base() >= rhs.as_base();
    }

#    ifdef __cpp_lib_three_way_comparison
    /*!\brief Performs a three-way comparison between `lhs` and `rhs`
     * \tparam UTypes The types of the elements of `rhs`. Automatically deduced.
     * \param lhs A tuple.
     * \param rhs A tuple with possibly different element types.
     * \returns An [ordering](https://en.cppreference.com/w/cpp/language/operator_comparison#Three-way_comparison)
     *          indicating the result of the comparison.
     */
    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes))
              && requires { requires (std::three_way_comparable_with<Types, UTypes> && ...); } // Defer instantiation
    constexpr friend auto operator<=>(tuple const & lhs, tuple<UTypes...> const & rhs)
    {
        return lhs.as_base() <=> rhs.as_base();
    }
#    endif // __cpp_lib_three_way_comparison
    //!\}

    /*!\name Comparison operators (std::tuple)
     * \{
     */
    /*!\brief Checks whether `lhs` and `rhs` are equal.
     * \tparam UTypes The types of the elements of `rhs`. Automatically deduced.
     * \param lhs A tuple.
     * \param rhs A std::tuple with possibly different element types.
     * \returns A bool indicating the result of the comparison.
     */
    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes))
              && requires { requires (std::equality_comparable_with<Types, UTypes> && ...); } // Defer instantiation
    constexpr friend bool operator==(tuple const & lhs, std::tuple<UTypes...> const & rhs)
    {
        return lhs.as_base() == rhs;
    }

    /*!\brief Checks whether `lhs` and `rhs` are unequal.
     * \tparam UTypes The types of the elements of `rhs`. Automatically deduced.
     * \param lhs A tuple.
     * \param rhs A std::tuple with possibly different element types.
     * \returns A bool indicating the result of the comparison.
     */
    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes))
              && requires { requires (std::equality_comparable_with<Types, UTypes> && ...); } // Defer instantiation
    constexpr friend bool operator!=(tuple const & lhs, std::tuple<UTypes...> const & rhs)
    {
        return lhs.as_base() != rhs;
    }

    /*!\brief Checks whether `lhs` is less than `rhs`.
     * \tparam UTypes The types of the elements of `rhs`. Automatically deduced.
     * \param lhs A tuple.
     * \param rhs A std::tuple with possibly different element types.
     * \returns A bool indicating the result of the comparison.
     */
    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes))
              && requires { requires (std::totally_ordered_with<Types, UTypes> && ...); } // Defer instantiation
    constexpr friend bool operator<(tuple const & lhs, std::tuple<UTypes...> const & rhs)
    {
        return lhs.as_base() < rhs;
    }

    /*!\brief Checks whether `lhs` is less than or equal to `rhs`.
     * \tparam UTypes The types of the elements of `rhs`. Automatically deduced.
     * \param lhs A tuple.
     * \param rhs A std::tuple with possibly different element types.
     * \returns A bool indicating the result of the comparison.
     */
    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes))
              && requires { requires (std::totally_ordered_with<Types, UTypes> && ...); } // Defer instantiation
    constexpr friend bool operator<=(tuple const & lhs, std::tuple<UTypes...> const & rhs)
    {
        return lhs.as_base() <= rhs;
    }

    /*!\brief Checks whether `lhs` is greater than `rhs`.
     * \tparam UTypes The types of the elements of `rhs`. Automatically deduced.
     * \param lhs A tuple.
     * \param rhs A std::tuple with possibly different element types.
     * \returns A bool indicating the result of the comparison.
     */
    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes))
              && requires { requires (std::totally_ordered_with<Types, UTypes> && ...); } // Defer instantiation
    constexpr friend bool operator>(tuple const & lhs, std::tuple<UTypes...> const & rhs)
    {
        return lhs.as_base() > rhs;
    }

    /*!\brief Checks whether `lhs` is greater than or equal to `rhs`.
     * \tparam UTypes The types of the elements of `rhs`. Automatically deduced.
     * \param lhs A tuple.
     * \param rhs A std::tuple with possibly different element types.
     * \returns A bool indicating the result of the comparison.
     */
    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes))
              && requires { requires (std::totally_ordered_with<Types, UTypes> && ...); } // Defer instantiation
    constexpr friend bool operator>=(tuple const & lhs, std::tuple<UTypes...> const & rhs)
    {
        return lhs.as_base() >= rhs;
    }

#    ifdef __cpp_lib_three_way_comparison
    /*!\brief Performs a three-way comparison between `lhs` and `rhs`
     * \tparam UTypes The types of the elements of `rhs`. Automatically deduced.
     * \param lhs A tuple.
     * \param rhs A std::tuple with possibly different element types.
     * \returns An [ordering](https://en.cppreference.com/w/cpp/language/operator_comparison#Three-way_comparison)
     *          indicating the result of the comparison.
     */
    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes))
              && requires { requires (std::three_way_comparable_with<Types, UTypes> && ...); } // Defer instantiation
    constexpr friend auto operator<=>(tuple const & lhs, std::tuple<UTypes...> const & rhs)
    {
        return lhs.as_base() <=> rhs;
    }
#    endif // __cpp_lib_three_way_comparison
    //!\}
};

//!\brief Class template argument deduction guide.
template <class... UTypes>
tuple(UTypes...) -> tuple<UTypes...>;

} // namespace seqan::stl

namespace seqan::stl::detail::tuple
{

template <typename query_t, typename... pack_t>
inline constexpr size_t count_in_pack = (std::is_same_v<query_t, pack_t> + ... + 0);

}

namespace std
{

template <typename... args>
struct tuple_size<seqan::stl::tuple<args...>> : public tuple_size<std::tuple<args...>>
{};

template <size_t index, typename... args>
struct tuple_element<index, seqan::stl::tuple<args...>> : public tuple_element<index, std::tuple<args...>>
{};

template <class... Ts, class... Us>
    requires requires { typename seqan::stl::tuple<std::common_type_t<Ts, Us>...>; }
struct common_type<seqan::stl::tuple<Ts...>, seqan::stl::tuple<Us...>>
{
    using type = seqan::stl::tuple<std::common_type_t<Ts, Us>...>;
};

template <class... Ts, class... Us>
    requires requires { typename seqan::stl::tuple<std::common_type_t<Ts, Us>...>; }
struct common_type<std::tuple<Ts...>, seqan::stl::tuple<Us...>>
{
    using type = seqan::stl::tuple<std::common_type_t<Ts, Us>...>;
};

template <class... Ts, class... Us>
    requires requires { typename seqan::stl::tuple<std::common_type_t<Ts, Us>...>; }
struct common_type<seqan::stl::tuple<Ts...>, std::tuple<Us...>>
{
    using type = seqan::stl::tuple<std::common_type_t<Ts, Us>...>;
};

template <class... Ts, class... Us, template <class> class TQual, template <class> class UQual>
    requires requires { typename seqan::stl::tuple<std::common_reference_t<TQual<Ts>, UQual<Us>>...>; }
struct basic_common_reference<seqan::stl::tuple<Ts...>, seqan::stl::tuple<Us...>, TQual, UQual>
{
    using type = seqan::stl::tuple<std::common_reference_t<TQual<Ts>, UQual<Us>>...>;
};

template <class... Ts, class... Us, template <class> class TQual, template <class> class UQual>
    requires requires { typename seqan::stl::tuple<std::common_reference_t<TQual<Ts>, UQual<Us>>...>; }
struct basic_common_reference<seqan::stl::tuple<Ts...>, std::tuple<Us...>, TQual, UQual>
{
    using type = seqan::stl::tuple<std::common_reference_t<TQual<Ts>, UQual<Us>>...>;
};

template <class... Ts, class... Us, template <class> class TQual, template <class> class UQual>
    requires requires { typename seqan::stl::tuple<std::common_reference_t<TQual<Ts>, UQual<Us>>...>; }
struct basic_common_reference<std::tuple<Ts...>, seqan::stl::tuple<Us...>, TQual, UQual>
{
    using type = seqan::stl::tuple<std::common_reference_t<TQual<Ts>, UQual<Us>>...>;
};

template <std::size_t i, typename... types>
constexpr std::tuple_element_t<i, seqan::stl::tuple<types...>> & get(seqan::stl::tuple<types...> & t) noexcept
    requires (i < sizeof...(types))
{
    return std::get<i>(static_cast<std::tuple<types...> &>(t));
}

template <std::size_t i, typename... types>
constexpr std::tuple_element_t<i, seqan::stl::tuple<types...>> const &
get(seqan::stl::tuple<types...> const & t) noexcept
    requires (i < sizeof...(types))
{
    return std::get<i>(static_cast<std::tuple<types...> const &>(t));
}

template <std::size_t i, typename... types>
constexpr std::tuple_element_t<i, seqan::stl::tuple<types...>> && get(seqan::stl::tuple<types...> && t) noexcept
    requires (i < sizeof...(types))
{
    return std::get<i>(static_cast<std::tuple<types...> &&>(std::move(t)));
}

template <std::size_t i, typename... types>
constexpr std::tuple_element_t<i, seqan::stl::tuple<types...>> const &&
get(seqan::stl::tuple<types...> const && t) noexcept
    requires (i < sizeof...(types))
{
    return std::get<i>(static_cast<std::tuple<types...> const &&>(std::move(t)));
}

template <typename type, typename... types>
constexpr type & get(seqan::stl::tuple<types...> & t) noexcept
    requires (seqan::stl::detail::tuple::count_in_pack<type, types...> == 1)
{
    return std::get<type>(static_cast<std::tuple<types...> &>(t));
}

template <typename type, typename... types>
constexpr type const & get(seqan::stl::tuple<types...> const & t) noexcept
    requires (seqan::stl::detail::tuple::count_in_pack<type, types...> == 1)
{
    return std::get<type>(static_cast<std::tuple<types...> const &>(t));
}

template <typename type, typename... types>
constexpr type && get(seqan::stl::tuple<types...> && t) noexcept
    requires (seqan::stl::detail::tuple::count_in_pack<type, types...> == 1)
{
    return std::get<type>(static_cast<std::tuple<types...> &&>(std::move(t)));
}

template <typename type, typename... types>
constexpr type const && get(seqan::stl::tuple<types...> const && t) noexcept
    requires (seqan::stl::detail::tuple::count_in_pack<type, types...> == 1)
{
    return std::get<type>(static_cast<std::tuple<types...> const &&>(std::move(t)));
}

} // namespace std

#endif // ifdef __cpp_lib_tuple_like

#endif // SEQAN_STD_TUPLE

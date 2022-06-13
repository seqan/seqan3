// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::common_tuple.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <tuple>

#include <seqan3/utility/tuple/common_pair.hpp>
#include <seqan3/utility/type_pack/traits.hpp>

namespace seqan3
{

/*!\brief A [std::tuple](https://en.cppreference.com/w/cpp/utility/tuple) implementation that incorporates most changes
 *        from C++23's standard library.
 * \ingroup utility_tuple
 */
template <class... Types>
class common_tuple : public std::tuple<Types...>
{
private:
    //!\brief The underlying std::tuple type.
    using base_t = std::tuple<Types...>;

    //!\brief Constructs by unfolding another tuple-like object.
    template <typename tuple_like_t, std::size_t... N>
    common_tuple(tuple_like_t && other, std::integer_sequence<size_t, N...>) :
        base_t((std::forward<std::tuple_element_t<N, tuple_like_t>>(std::get<N>(other)))...)
    {}

public:
    /*!\name Default constructors and assignments.
     * \{
     */
    common_tuple() = default;                                 //!< Defaulted.
    common_tuple(common_tuple const &) = default;             //!< Defaulted.
    common_tuple & operator=(common_tuple const &) = default; //!< Defaulted.
    ~common_tuple() = default;                                //!< Defaulted.
    //!\}

    /*!\brief Returns the common_tuple as the underlying std::tuple type.
     * \private
     */
    base_t & as_base() noexcept
    {
        return *this;
    }

    /*!\brief Returns the common_tuple as the underlying std::tuple type.
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
    constexpr common_tuple(UTypes &... other) : base_t(other...)
    {}

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_constructible_v<Types, UTypes const> && ...)
    constexpr common_tuple(UTypes const &... other) : base_t(other...)
    {}

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_constructible_v<Types, UTypes> && ...)
    constexpr common_tuple(UTypes &&... other) : base_t(std::forward<UTypes>(other)...)
    {}

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_constructible_v<Types, UTypes const> && ...)
    constexpr common_tuple(UTypes const &&... other) : base_t(std::forward<UTypes const>(other)...)
    {}
    //!@}
    //!\}

    /*!\name Construct from common_tuple.
     * \{
     */
    //!@{ Constructs from common_tuple.
    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_constructible_v<Types, UTypes &> && ...)
    constexpr common_tuple(common_tuple<UTypes...> & other) : common_tuple(other, std::index_sequence_for<Types...>{})
    {}

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_constructible_v<Types, UTypes const> && ...)
    constexpr common_tuple(common_tuple<UTypes...> const & other) :
        common_tuple(other, std::index_sequence_for<Types...>{})
    {}

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_constructible_v<Types, UTypes> && ...)
    constexpr common_tuple(common_tuple<UTypes...> && other) :
        common_tuple(std::move(other), std::index_sequence_for<Types...>{})
    {}

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_constructible_v<Types, UTypes const> && ...)
    constexpr common_tuple(common_tuple<UTypes...> const && other) :
        common_tuple(std::move(other), std::index_sequence_for<Types...>{})
    {}
    //!@}
    //!\}

    /*!\name Construct from common_pair.
     * \{
     */
    //!@{ Constructs from common_pair.
    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_constructible_v<Types, UTypes &> && ...)
    constexpr common_tuple(common_pair<UTypes...> & other) : common_tuple(other, std::index_sequence_for<Types...>{})
    {}

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_constructible_v<Types, UTypes const> && ...)
    constexpr common_tuple(common_pair<UTypes...> const & other) :
        common_tuple(other, std::index_sequence_for<Types...>{})
    {}

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_constructible_v<Types, UTypes> && ...)
    constexpr common_tuple(common_pair<UTypes...> && other) :
        common_tuple(std::move(other), std::index_sequence_for<Types...>{})
    {}

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_constructible_v<Types, UTypes const> && ...)
    constexpr common_tuple(common_pair<UTypes...> const && other) :
        common_tuple(std::move(other), std::index_sequence_for<Types...>{})
    {}
    //!@}
    //!\}

    /*!\name Construct from std::tuple.
     * \{
     */
    //!@{ Constructs from std::tuple.
    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_constructible_v<Types, UTypes &> && ...)
    constexpr common_tuple(std::tuple<UTypes...> & other) : common_tuple(other, std::index_sequence_for<Types...>{})
    {}

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_constructible_v<Types, UTypes const> && ...)
    constexpr common_tuple(std::tuple<UTypes...> const & other) :
        common_tuple(other, std::index_sequence_for<Types...>{})
    {}

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_constructible_v<Types, UTypes> && ...)
    constexpr common_tuple(std::tuple<UTypes...> && other) :
        common_tuple(std::move(other), std::index_sequence_for<Types...>{})
    {}

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_constructible_v<Types, UTypes const> && ...)
    constexpr common_tuple(std::tuple<UTypes...> const && other) :
        common_tuple(std::move(other), std::index_sequence_for<Types...>{})
    {}

    /*!\name Construct from std::pair.
     * \{
     */
    //!@{ Constructs from std::pair.
    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_constructible_v<Types, UTypes &> && ...)
    constexpr common_tuple(std::pair<UTypes...> & other) : common_tuple(other, std::index_sequence_for<Types...>{})
    {}

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_constructible_v<Types, UTypes const> && ...)
    constexpr common_tuple(std::pair<UTypes...> const & other) :
        common_tuple(other, std::index_sequence_for<Types...>{})
    {}

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_constructible_v<Types, UTypes> && ...)
    constexpr common_tuple(std::pair<UTypes...> && other) :
        common_tuple(std::move(other), std::index_sequence_for<Types...>{})
    {}

    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_constructible_v<Types, UTypes const> && ...)
    constexpr common_tuple(std::pair<UTypes...> const && other) :
        common_tuple(std::move(other), std::index_sequence_for<Types...>{})
    {}
    //!@}
    //!\}

    /*!\name Assign from common_tuple.
     * \{
     */
    //!@{ Assigns from common_tuple.
    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_assignable_v<Types, UTypes &> && ...)
    constexpr common_tuple & operator=(common_tuple<UTypes...> & other)
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
    constexpr common_tuple & operator=(common_tuple<UTypes...> const & other)
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
    constexpr common_tuple & operator=(common_tuple<UTypes...> && other)
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
    constexpr common_tuple & operator=(common_tuple<UTypes...> const && other)
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
    constexpr common_tuple const & operator=(common_tuple<UTypes...> & other) const
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
    constexpr common_tuple const & operator=(common_tuple<UTypes...> const & other) const
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
    constexpr common_tuple const & operator=(common_tuple<UTypes...> && other) const
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
    constexpr common_tuple const & operator=(common_tuple<UTypes...> const && other) const
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

    /*!\name Assign from common_pair.
     * \{
     */
    //!@{ Assigns from common_pair.
    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes)) && (std::is_assignable_v<Types, UTypes &> && ...)
    constexpr common_tuple & operator=(common_pair<UTypes...> & other)
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
    constexpr common_tuple & operator=(common_pair<UTypes...> const & other)
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
    constexpr common_tuple & operator=(common_pair<UTypes...> && other)
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
    constexpr common_tuple & operator=(common_pair<UTypes...> const && other)
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
    constexpr common_tuple const & operator=(common_pair<UTypes...> & other) const
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
    constexpr common_tuple const & operator=(common_pair<UTypes...> const & other) const
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
    constexpr common_tuple const & operator=(common_pair<UTypes...> && other) const
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
    constexpr common_tuple const & operator=(common_pair<UTypes...> const && other) const
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
    constexpr common_tuple & operator=(std::tuple<UTypes...> & other)
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
    constexpr common_tuple & operator=(std::tuple<UTypes...> const & other)
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
    constexpr common_tuple & operator=(std::tuple<UTypes...> && other)
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
    constexpr common_tuple & operator=(std::tuple<UTypes...> const && other)
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
    constexpr common_tuple const & operator=(std::tuple<UTypes...> & other) const
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
    constexpr common_tuple const & operator=(std::tuple<UTypes...> const & other) const
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
    constexpr common_tuple const & operator=(std::tuple<UTypes...> && other) const
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
    constexpr common_tuple const & operator=(std::tuple<UTypes...> const && other) const
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
    constexpr common_tuple & operator=(std::pair<UTypes...> & other)
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
    constexpr common_tuple & operator=(std::pair<UTypes...> const & other)
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
    constexpr common_tuple & operator=(std::pair<UTypes...> && other)
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
    constexpr common_tuple & operator=(std::pair<UTypes...> const && other)
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
    constexpr common_tuple const & operator=(std::pair<UTypes...> & other) const
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
    constexpr common_tuple const & operator=(std::pair<UTypes...> const & other) const
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
    constexpr common_tuple const & operator=(std::pair<UTypes...> && other) const
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
    constexpr common_tuple const & operator=(std::pair<UTypes...> const && other) const
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

    /*!\name Comparison operators (common_tuple)
     * \{
     */
    /*!\brief Checks whether `lhs` and `rhs` are equal.
     * \tparam UTypes The types of the elements of `rhs`. Automatically deduced.
     * \param lhs A common_tuple.
     * \param rhs A common_tuple with possibly different element types.
     * \returns A bool indicating the result of the comparison.
     */
    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes))
              && requires { requires (std::equality_comparable_with<Types, UTypes> && ...); } // Defer instantiation
    constexpr friend bool operator==(common_tuple const & lhs, common_tuple<UTypes...> const & rhs)
    {
        static_assert((std::equality_comparable_with<Types, UTypes> && ...));
        return lhs.as_base() == rhs.as_base();
    }

    /*!\brief Checks whether `lhs` and `rhs` are unequal.
     * \tparam UTypes The types of the elements of `rhs`. Automatically deduced.
     * \param lhs A common_tuple.
     * \param rhs A common_tuple with possibly different element types.
     * \returns A bool indicating the result of the comparison.
     */
    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes))
              && requires { requires (std::equality_comparable_with<Types, UTypes> && ...); } // Defer instantiation
    constexpr friend bool operator!=(common_tuple const & lhs, common_tuple<UTypes...> const & rhs)
    {
        return lhs.as_base() != rhs.as_base();
    }

    /*!\brief Checks whether `lhs` is less than `rhs`.
     * \tparam UTypes The types of the elements of `rhs`. Automatically deduced.
     * \param lhs A common_tuple.
     * \param rhs A common_tuple with possibly different element types.
     * \returns A bool indicating the result of the comparison.
     */
    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes))
              && requires { requires (std::totally_ordered_with<Types, UTypes> && ...); } // Defer instantiation
    constexpr friend bool operator<(common_tuple const & lhs, common_tuple<UTypes...> const & rhs)
    {
        return lhs.as_base() < rhs.as_base();
    }

    /*!\brief Checks whether `lhs` is less than or equal to `rhs`.
     * \tparam UTypes The types of the elements of `rhs`. Automatically deduced.
     * \param lhs A common_tuple.
     * \param rhs A common_tuple with possibly different element types.
     * \returns A bool indicating the result of the comparison.
     */
    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes))
              && requires { requires (std::totally_ordered_with<Types, UTypes> && ...); } // Defer instantiation
    constexpr friend bool operator<=(common_tuple const & lhs, common_tuple<UTypes...> const & rhs)
    {
        return lhs.as_base() <= rhs.as_base();
    }

    /*!\brief Checks whether `lhs` is greater than `rhs`.
     * \tparam UTypes The types of the elements of `rhs`. Automatically deduced.
     * \param lhs A common_tuple.
     * \param rhs A common_tuple with possibly different element types.
     * \returns A bool indicating the result of the comparison.
     */
    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes))
              && requires { requires (std::totally_ordered_with<Types, UTypes> && ...); } // Defer instantiation
    constexpr friend bool operator>(common_tuple const & lhs, common_tuple<UTypes...> const & rhs)
    {
        return lhs.as_base() > rhs.as_base();
    }

    /*!\brief Checks whether `lhs` is greater than or equal to `rhs`.
     * \tparam UTypes The types of the elements of `rhs`. Automatically deduced.
     * \param lhs A common_tuple.
     * \param rhs A common_tuple with possibly different element types.
     * \returns A bool indicating the result of the comparison.
     */
    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes))
              && requires { requires (std::totally_ordered_with<Types, UTypes> && ...); } // Defer instantiation
    constexpr friend bool operator>=(common_tuple const & lhs, common_tuple<UTypes...> const & rhs)
    {
        return lhs.as_base() >= rhs.as_base();
    }

#ifdef __cpp_lib_three_way_comparison
    /*!\brief Performs a three-way comparison between `lhs` and `rhs`
     * \tparam UTypes The types of the elements of `rhs`. Automatically deduced.
     * \param lhs A common_tuple.
     * \param rhs A common_tuple with possibly different element types.
     * \returns An [ordering](https://en.cppreference.com/w/cpp/language/operator_comparison#Three-way_comparison)
     *          indicating the result of the comparison.
     */
    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes))
              && requires { requires (std::three_way_comparable_with<Types, UTypes> && ...); } // Defer instantiation
    constexpr friend auto operator<=>(common_tuple const & lhs, common_tuple<UTypes...> const & rhs)
    {
        return lhs.as_base() <=> rhs.as_base();
    }
#endif // __cpp_lib_three_way_comparison
    //!\}

    /*!\name Comparison operators (std::tuple)
     * \{
     */
    /*!\brief Checks whether `lhs` and `rhs` are equal.
     * \tparam UTypes The types of the elements of `rhs`. Automatically deduced.
     * \param lhs A common_tuple.
     * \param rhs A std::tuple with possibly different element types.
     * \returns A bool indicating the result of the comparison.
     */
    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes))
              && requires { requires (std::equality_comparable_with<Types, UTypes> && ...); } // Defer instantiation
    constexpr friend bool operator==(common_tuple const & lhs, std::tuple<UTypes...> const & rhs)
    {
        return lhs.as_base() == rhs;
    }

    /*!\brief Checks whether `lhs` and `rhs` are unequal.
     * \tparam UTypes The types of the elements of `rhs`. Automatically deduced.
     * \param lhs A common_tuple.
     * \param rhs A std::tuple with possibly different element types.
     * \returns A bool indicating the result of the comparison.
     */
    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes))
              && requires { requires (std::equality_comparable_with<Types, UTypes> && ...); } // Defer instantiation
    constexpr friend bool operator!=(common_tuple const & lhs, std::tuple<UTypes...> const & rhs)
    {
        return lhs.as_base() != rhs;
    }

    /*!\brief Checks whether `lhs` is less than `rhs`.
     * \tparam UTypes The types of the elements of `rhs`. Automatically deduced.
     * \param lhs A common_tuple.
     * \param rhs A std::tuple with possibly different element types.
     * \returns A bool indicating the result of the comparison.
     */
    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes))
              && requires { requires (std::totally_ordered_with<Types, UTypes> && ...); } // Defer instantiation
    constexpr friend bool operator<(common_tuple const & lhs, std::tuple<UTypes...> const & rhs)
    {
        return lhs.as_base() < rhs;
    }

    /*!\brief Checks whether `lhs` is less than or equal to `rhs`.
     * \tparam UTypes The types of the elements of `rhs`. Automatically deduced.
     * \param lhs A common_tuple.
     * \param rhs A std::tuple with possibly different element types.
     * \returns A bool indicating the result of the comparison.
     */
    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes))
              && requires { requires (std::totally_ordered_with<Types, UTypes> && ...); } // Defer instantiation
    constexpr friend bool operator<=(common_tuple const & lhs, std::tuple<UTypes...> const & rhs)
    {
        return lhs.as_base() <= rhs;
    }

    /*!\brief Checks whether `lhs` is greater than `rhs`.
     * \tparam UTypes The types of the elements of `rhs`. Automatically deduced.
     * \param lhs A common_tuple.
     * \param rhs A std::tuple with possibly different element types.
     * \returns A bool indicating the result of the comparison.
     */
    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes))
              && requires { requires (std::totally_ordered_with<Types, UTypes> && ...); } // Defer instantiation
    constexpr friend bool operator>(common_tuple const & lhs, std::tuple<UTypes...> const & rhs)
    {
        return lhs.as_base() > rhs;
    }

    /*!\brief Checks whether `lhs` is greater than or equal to `rhs`.
     * \tparam UTypes The types of the elements of `rhs`. Automatically deduced.
     * \param lhs A common_tuple.
     * \param rhs A std::tuple with possibly different element types.
     * \returns A bool indicating the result of the comparison.
     */
    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes))
              && requires { requires (std::totally_ordered_with<Types, UTypes> && ...); } // Defer instantiation
    constexpr friend bool operator>=(common_tuple const & lhs, std::tuple<UTypes...> const & rhs)
    {
        return lhs.as_base() >= rhs;
    }

#ifdef __cpp_lib_three_way_comparison
    /*!\brief Performs a three-way comparison between `lhs` and `rhs`
     * \tparam UTypes The types of the elements of `rhs`. Automatically deduced.
     * \param lhs A common_tuple.
     * \param rhs A std::tuple with possibly different element types.
     * \returns An [ordering](https://en.cppreference.com/w/cpp/language/operator_comparison#Three-way_comparison)
     *          indicating the result of the comparison.
     */
    template <class... UTypes>
        requires (sizeof...(Types) == sizeof...(UTypes))
              && requires { requires (std::three_way_comparable_with<Types, UTypes> && ...); } // Defer instantiation
    constexpr friend auto operator<=>(common_tuple const & lhs, std::tuple<UTypes...> const & rhs)
    {
        return lhs.as_base() <=> rhs;
    }
#endif // __cpp_lib_three_way_comparison
    //!\}
};

//!\brief Class template argument deduction guide.
template <class... UTypes>
common_tuple(UTypes...) -> common_tuple<UTypes...>;

} // namespace seqan3

//!\cond
namespace std
{

template <typename... args>
struct tuple_size<seqan3::common_tuple<args...>> : public tuple_size<std::tuple<args...>>
{};

template <size_t index, typename... args>
struct tuple_element<index, seqan3::common_tuple<args...>> : public tuple_element<index, std::tuple<args...>>
{};

template <class... Ts, class... Us>
    requires requires { typename seqan3::common_tuple<std::common_type_t<Ts, Us>...>; }
struct common_type<seqan3::common_tuple<Ts...>, seqan3::common_tuple<Us...>>
{
    using type = seqan3::common_tuple<std::common_type_t<Ts, Us>...>;
};

template <class... Ts, class... Us>
    requires requires { typename seqan3::common_tuple<std::common_type_t<Ts, Us>...>; }
struct common_type<std::tuple<Ts...>, seqan3::common_tuple<Us...>>
{
    using type = seqan3::common_tuple<std::common_type_t<Ts, Us>...>;
};

template <class... Ts, class... Us>
    requires requires { typename seqan3::common_tuple<std::common_type_t<Ts, Us>...>; }
struct common_type<seqan3::common_tuple<Ts...>, std::tuple<Us...>>
{
    using type = seqan3::common_tuple<std::common_type_t<Ts, Us>...>;
};

template <class... Ts, class... Us, template <class> class TQual, template <class> class UQual>
    requires requires { typename seqan3::common_tuple<std::common_reference_t<TQual<Ts>, UQual<Us>>...>; }
struct basic_common_reference<seqan3::common_tuple<Ts...>, seqan3::common_tuple<Us...>, TQual, UQual>
{
    using type = seqan3::common_tuple<std::common_reference_t<TQual<Ts>, UQual<Us>>...>;
};

template <class... Ts, class... Us, template <class> class TQual, template <class> class UQual>
    requires requires { typename seqan3::common_tuple<std::common_reference_t<TQual<Ts>, UQual<Us>>...>; }
struct basic_common_reference<seqan3::common_tuple<Ts...>, std::tuple<Us...>, TQual, UQual>
{
    using type = seqan3::common_tuple<std::common_reference_t<TQual<Ts>, UQual<Us>>...>;
};

template <class... Ts, class... Us, template <class> class TQual, template <class> class UQual>
    requires requires { typename seqan3::common_tuple<std::common_reference_t<TQual<Ts>, UQual<Us>>...>; }
struct basic_common_reference<std::tuple<Ts...>, seqan3::common_tuple<Us...>, TQual, UQual>
{
    using type = seqan3::common_tuple<std::common_reference_t<TQual<Ts>, UQual<Us>>...>;
};

template <std::size_t i, typename... types>
constexpr std::tuple_element_t<i, seqan3::common_tuple<types...>> & get(seqan3::common_tuple<types...> & t) noexcept
    requires (i < sizeof...(types))
{
    return std::get<i>(static_cast<std::tuple<types...> &>(t));
}

template <std::size_t i, typename... types>
constexpr std::tuple_element_t<i, seqan3::common_tuple<types...>> const &
get(seqan3::common_tuple<types...> const & t) noexcept
    requires (i < sizeof...(types))
{
    return std::get<i>(static_cast<std::tuple<types...> const &>(t));
}

template <std::size_t i, typename... types>
constexpr std::tuple_element_t<i, seqan3::common_tuple<types...>> && get(seqan3::common_tuple<types...> && t) noexcept
    requires (i < sizeof...(types))
{
    return std::get<i>(static_cast<std::tuple<types...> &&>(std::move(t)));
}

template <std::size_t i, typename... types>
constexpr std::tuple_element_t<i, seqan3::common_tuple<types...>> const &&
get(seqan3::common_tuple<types...> const && t) noexcept
    requires (i < sizeof...(types))
{
    return std::get<i>(static_cast<std::tuple<types...> const &&>(std::move(t)));
}

template <typename type, typename... types>
constexpr type & get(seqan3::common_tuple<types...> & t) noexcept
    requires (seqan3::pack_traits::count<type, types...> == 1)
{
    return std::get<type>(static_cast<std::tuple<types...> &>(t));
}

template <typename type, typename... types>
constexpr type const & get(seqan3::common_tuple<types...> const & t) noexcept
    requires (seqan3::pack_traits::count<type, types...> == 1)
{
    return std::get<type>(static_cast<std::tuple<types...> const &>(t));
}

template <typename type, typename... types>
constexpr type && get(seqan3::common_tuple<types...> && t) noexcept
    requires (seqan3::pack_traits::count<type, types...> == 1)
{
    return std::get<type>(static_cast<std::tuple<types...> &&>(std::move(t)));
}

template <typename type, typename... types>
constexpr type const && get(seqan3::common_tuple<types...> const && t) noexcept
    requires (seqan3::pack_traits::count<type, types...> == 1)
{
    return std::get<type>(static_cast<std::tuple<types...> const &&>(std::move(t)));
}

} // namespace std
//!\endcond

// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
//
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================

/*!\file
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 * \author David Heller <david.heller AT fu-berlin.de>
 * \brief Contains seqan3::union_composition.
 */

#pragma once

#include <algorithm>
#include <array>
#include <utility>
#include <cassert>
#include <variant>

#include <meta/meta.hpp>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/core/detail/int_types.hpp>
#include <seqan3/std/concepts>

namespace seqan3
{

/*!\brief A combined alphabet that can hold values of either of its alternatives.
 * \ingroup composition
 * \tparam ...alternative_types Types of possible values (at least 2); all must model seqan3::alphabet_concept and be
 *                              unique.
 * \implements seqan3::alphabet_concept
 * \implements seqan3::detail::constexpr_alphabet_concept
 *
 * \details
 *
 * The union_composition represents the union of two or more alternative alphabets (e.g. the
 * four letter DNA alternative + the gap alternative). It behaves similar to a
 * [union](https://en.cppreference.com/w/cpp/language/union) or std::variant, but it preserves the
 * seqan3::alphabet_concept.
 *
 * Short description:
 *   * combines multiple different alphabets in an "either-or"-fashion;
 *   * is itself a seqan3::alphabet_concept;
 *   * its alphabet size is the sum of the individual sizes;
 *   * default initialises to the the first alternative's default (no empty state like std::variant);
 *   * constructible, assignable and (in-)equality-comparable with each alternative type and also all types that
 *     these are constructible/assignable/equality-comparable with;
 *   * only convertible to its alternatives through the member function convert_to() (which can throw!)
 *
 * ### Example
 *
 * \snippet test/snippet/alphabet/composition/union_composition.cpp usage
 */
template <typename ...alternative_types>
//!\cond
    requires (detail::constexpr_alphabet_concept<alternative_types> && ... && (sizeof...(alternative_types) >= 2))
//!\endcond
class union_composition
{
private:
    //!\brief A meta::list of the types of each alternative in the composition
    using alternatives = meta::list<alternative_types...>;

    static_assert(std::Same<alternatives, meta::unique<alternatives>>,
                  "All types in a union_composition must be distinct.");

    /*!\brief 'Callable' helper class that is invokable by meta::invoke.
      * Returns an std::true_type if the `type` is constructable from `T`.
      */
    template <typename T>
    struct constructible_from
    {
        //!\brief The returned type when invoked.
        template <typename type>
        using invoke = std::integral_constant<bool, std::is_constructible_v<type, T>>;
    };

    /*!\brief 'Callable' helper class that is invokable by meta::invoke.
     * Returns an std::true_type if the `T` is implicitly convertible to `type`.
     */
    template <typename T>
    struct implicitly_convertible_from
    {
        //!\brief The returned type when invoked.
        template <typename type>
        using invoke = std::integral_constant<bool, implicitly_convertible_to_concept<T, type>>;
    };

    /*!\brief 'Callable' helper class that is invokable by meta::invoke.
     * Returns an std::true_type if the `type` is assignable from `T`.
     */
    template <typename T>
    struct assignable_from
    {
        //!\brief The returned type when invoked.
        template <typename type>
        using invoke = std::integral_constant<bool, std::Assignable<type, T>>;
    };

    /*!\brief 'Callable' helper class that is invokable by meta::invoke.
     * Returns an std::true_type if the `type` is weakly equality comparable to `T`.
     */
    template <typename T>
    struct weakly_equality_comparable_with
    {
        //!\brief The returned type when invoked.
        template <typename type>
        using invoke = std::integral_constant<bool, std::detail::WeaklyEqualityComparableWith<type, T>>;
    };

    //!\brief Ridiculously verbose helper, because we don't have local concept definitions or template specialisations.
    template <template <typename> typename fun_t, typename indirect_alternative_t>
    constexpr static bool one_alternative_is_helper()
    {
        // filter out this and "direct" alternatives
        if constexpr (std::is_same_v<indirect_alternative_t, union_composition> ||
                      meta::in<alternatives, indirect_alternative_t>::value)
        {
            return false;
        }
        else // make the second term only be instantiated if the prior condition is not met to prevent recursions
        {
            return !meta::empty<meta::find_if<alternatives, fun_t<indirect_alternative_t>>>::value;
        }
    }

    /*!\brief 'Is set to `true` if one alternative type in `alternatives` evaluates
     * to true when invoked with `indirect_alternative_t` by fun_t.
     */
    template <template <typename> typename fun_t, typename indirect_alternative_t>
    static constexpr bool one_alternative_is = one_alternative_is_helper<fun_t, indirect_alternative_t>();

public:
    /*!\brief Returns true if alternative_t is one of the given alternative types.
     * \tparam alternative_t The type to check.
     *
     * \snippet test/snippet/alphabet/composition/union_composition.cpp holds_alternative
     */
    template <typename alternative_t>
    static constexpr bool holds_alternative() noexcept
    {
        return meta::in<alternatives, alternative_t>::value;
    }

    //!\brief The size of the alternative, i.e. the number of different values it can take.
    static constexpr auto value_size =
        detail::min_viable_uint_v<(static_cast<size_t>(alphabet_size_v<alternative_types>) + ...)>;

    //!\brief The type of the alternative when converted to char (e.g. via \link to_char \endlink).
    using char_type = underlying_char_t<meta::front<alternatives>>;

    //!\brief The type of the alternative when represented as a number (e.g. via \link to_rank \endlink).
    using rank_type = detail::min_viable_uint_t<value_size>;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr union_composition() = default;
    constexpr union_composition(union_composition const &) = default;
    constexpr union_composition(union_composition &&) = default;
    constexpr union_composition & operator=(union_composition const &) = default;
    constexpr union_composition & operator=(union_composition &&) = default;

    /*!\brief Construction via the value of an alternative.
     * \tparam alternative_t One of the alternative types.
     * \param  alternative   The value of a alternative that should be assigned.
     *
     * \snippet test/snippet/alphabet/composition/union_composition.cpp value construction
     */
    template <typename alternative_t>
    //!\cond
        requires holds_alternative<alternative_t>()
    //!\endcond
    constexpr union_composition(alternative_t const & alternative) noexcept :
        _value{rank_by_type_(alternative)}
    {}

    /*!\brief Construction via the value of a type that an alternative type is constructible from.
     * \tparam indirect_alternative_t A type that one of the alternative types is constructible from.
     * \param  rhs The value that should be assigned.
     *
     * \snippet test/snippet/alphabet/composition/union_composition.cpp conversion
     * \attention When selecting the alternative alphabet types which require only implicit conversion
     * or constructor calls, are preferred over those that require explicit ones.
     */
    template <typename indirect_alternative_t>
    //!\cond
        requires !one_alternative_is<implicitly_convertible_from, indirect_alternative_t> &&
                 one_alternative_is<constructible_from, indirect_alternative_t>
    //!\endcond
    constexpr union_composition(indirect_alternative_t const & rhs) noexcept :
        _value{rank_by_type_(meta::front<meta::find_if<alternatives,
                                                       constructible_from<indirect_alternative_t>>>(rhs))}
    {}

    //!\cond
    template <typename indirect_alternative_t>
        requires one_alternative_is<implicitly_convertible_from, indirect_alternative_t>
    constexpr union_composition(indirect_alternative_t const & rhs) noexcept :
        _value{rank_by_type_(meta::front<meta::find_if<alternatives,
                                                       implicitly_convertible_from<indirect_alternative_t>>>(rhs))}
    {}
    //!\endcond

    /*!\brief Assignment via a value that one of the alternative types is assignable from.
     * \tparam indirect_alternative_t A type that one of the alternatives is assignable from.
     * \param  rhs The value of an alternative.
     *
     * \snippet test/snippet/alphabet/composition/union_composition.cpp subtype_construction
     */
    template <typename indirect_alternative_t>
    //!\cond
        requires one_alternative_is<assignable_from, indirect_alternative_t>
    //!\endcond
    constexpr union_composition & operator=(indirect_alternative_t const & rhs) noexcept
    {
        using alternative_t = meta::front<meta::find_if<alternatives, assignable_from<indirect_alternative_t>>>;
        alternative_t alternative = rhs;
        _value = rank_by_type_(alternative);
        return *this;
    }
    //!\}

    /*!\name Read functions
     * \{
     */
    //!\brief Return the letter as a character of char_type.
    constexpr char_type to_char() const noexcept
    {
        return value_to_char[_value];
    }

    //!\brief Return the letter's numeric value or rank in the alternative.
    constexpr rank_type to_rank() const noexcept
    {
        return _value;
    }
    //!\}

    /*!\name Write functions
     * \{
     */
    //!\brief Assign from a character.
    //!\param c The character to assign to.
    constexpr union_composition & assign_char(char_type const c) noexcept
    {
        using index_t = std::make_unsigned_t<char_type>;
        _value = char_to_value[static_cast<index_t>(c)];
        return *this;
    }

    //!\brief Assign from a numeric value.
    //!\param i The rank of a character to assign to.
    constexpr union_composition & assign_rank(rank_type const i) noexcept
    {
        assert(i < value_size);
        _value = i;
        return *this;
    }
    //!\}

    /*!\name Conversion (by index)
     * \{
     */
    //!\brief Whether the union alphabet currently holds a value of the given alternative.
    //!\tparam index Index of the alternative to check for.
    template <size_t index>
    constexpr bool is_alternative() const noexcept
    {
        static_assert(index < value_size, "The union_composition contains less alternatives than you are checking.");
        return (to_rank() >= partial_sum_sizes[index]) && (to_rank() < partial_sum_sizes[index + 1]);
    }

    /*!\brief Convert to the specified alphabet (throws if is_alternative() would be false).
     * \tparam index Index of the alternative to check for.
     * \throws std::bad_variant_access If the union_alphabet currently holds the value of a different alternative.
     */
    template <size_t index>
    constexpr auto convert_to() const
    {
        return convert_impl<index, true>();
    }

    /*!\brief Convert to the specified alphabet (**undefined behaviour** if is_alternative() would be false).
     * \tparam index Index of the alternative to check for.
     */
    template <size_t index>
    constexpr auto convert_unsafely_to() const noexcept
    {
        return convert_impl<index, false>();
    }
    //!\}

    /*!\name Conversion (by type)
     * \{
     */
    /*!\copybrief is_alternative()
     * \tparam alternative_t The type of the alternative that you wish to check for.
     */
    template <typename alternative_t>
    constexpr bool is_alternative() const noexcept
        requires holds_alternative<alternative_t>()
    {
        constexpr size_t index = meta::find_index<alternatives, alternative_t>::value;
        return is_alternative<index>();
    }

    /*!\copybrief convert_to()
     * \tparam alternative_t The type of the alternative that you wish to check for.
     * \throws std::bad_variant_access If the union_alphabet currently holds the value of a different alternative.
     */
    template <typename alternative_t>
    constexpr alternative_t convert_to() const
        requires holds_alternative<alternative_t>()
    {
        constexpr size_t index = meta::find_index<alternatives, alternative_t>::value;
        return convert_impl<index, true>();
    }

    /*!\copybrief convert_unsafely_to()
     * \tparam alternative_t The type of the alternative that you wish to check for.
     */
    template <typename alternative_t>
    constexpr alternative_t convert_unsafely_to() const noexcept
        requires holds_alternative<alternative_t>()
    {
        constexpr size_t index = meta::find_index<alternatives, alternative_t>::value;
        return convert_impl<index, false>();
    }
    //!\}

    /*!\name Comparison operators (against self)
     * \{
     */
    constexpr bool operator==(union_composition const & rhs) const noexcept
    {
        return _value == rhs._value;
    }

    constexpr bool operator!=(union_composition const & rhs) const noexcept
    {
        return _value != rhs._value;
    }

    constexpr bool operator<(union_composition const & rhs) const noexcept
    {
        return _value < rhs._value;
    }

    constexpr bool operator>(union_composition const & rhs) const noexcept
    {
        return _value > rhs._value;
    }

    constexpr bool operator<=(union_composition const & rhs) const noexcept
    {
        return _value <= rhs._value;
    }

    constexpr bool operator>=(union_composition const & rhs) const noexcept
    {
        return _value >= rhs._value;
    }
    //!\}

    /*!\name Comparison operators (against alternatives)
     * \brief Defines comparison against alternatives, e.g. `union_composition<dna5, gap>{gap::GAP} == dna5::C`. Only
     *        (in-)equality comparison is explicitly defined, because it would be difficult to argue about e.g.
     *        `union_composition<dna5, gap>{gap::GAP} < dna5::C`.
     * \{
     */
    template <typename alternative_t>
    constexpr bool operator==(alternative_t const & rhs) const noexcept
        requires holds_alternative<alternative_t>()
    {
        return is_alternative<alternative_t>() && (convert_unsafely_to<alternative_t>() == rhs);
    }

    template <typename alternative_t>
    constexpr bool operator!=(alternative_t const & rhs) const noexcept
        requires holds_alternative<alternative_t>()
    {
        return !operator==(rhs);
    }
    //!\}

    /*!\name Comparison operators (against indirect alternatives)
     * \brief Defines comparison against types that are comparable with alternatives, e.g.
     *        `union_composition<dna5, gap>{dna5::C} == rna5::C`. Only (in-)equality comparison is explicitly defined,
     *        because it would be difficult to argue about e.g.
     *        `union_composition<dna5, gap>{gap::GAP} < rna5::C`.
     * \{
     */
    template <typename indirect_alternative_type>
    constexpr bool operator==(indirect_alternative_type const & rhs) const noexcept
    //!\cond
        requires one_alternative_is<weakly_equality_comparable_with, indirect_alternative_type>
    //!\endcond
    {
        using alternative_t =
            meta::front<meta::find_if<alternatives, weakly_equality_comparable_with<indirect_alternative_type>>>;
        return is_alternative<alternative_t>() && (convert_unsafely_to<alternative_t>() == rhs);
    }

    template <typename indirect_alternative_type>
    constexpr bool operator!=(indirect_alternative_type const & rhs) const noexcept
    //!\cond
        requires one_alternative_is<weakly_equality_comparable_with, indirect_alternative_type>
    //!\endcond
    {
        return !operator==(rhs);
    }
    //!\}

    //!\privatesection
    //!\brief The data member.
    rank_type _value;
    //!\publicsection

protected:
    //!\privatesection

    /*!\brief Implementation function for convert_to() and convert_unsafely_to().
     * \tparam index  Index of the alternative to convert to.
     * \tparam throws Whether to perform checks (and throw) or not.
     */
    template <size_t index, bool throws>
    constexpr auto convert_impl() const noexcept(!throws) -> meta::at_c<alternatives, index>
    {
        static_assert(index < value_size, "The union_composition contains less alternatives than you are checking.");
        using alternative_t = meta::at_c<alternatives, index>;

        if constexpr (throws)
        {
            if (!is_alternative<index>()) // [[unlikely]]
            {
                throw std::bad_variant_access{};
            }
        }

        using seqan3::assign_rank;
        return assign_rank(alternative_t{}, to_rank() - partial_sum_sizes[index]);
    }

    /*!\brief Compile-time generated lookup table which contains the partial
     * sum up to the position of each alternative.
     *
     * An array which contains the prefix sum over all
     * alternative_types::value_size's.
     *
     */
    static constexpr std::array partial_sum_sizes = []() constexpr
    {
        constexpr size_t N = sizeof...(alternative_types) + 1;

        std::array<rank_type, N> partial_sum{0, alphabet_size_v<alternative_types>...};
        for (size_t i = 1u; i < N; ++i)
            partial_sum[i] += partial_sum[i-1];

        return partial_sum;
    }();

    /*!\brief Compile-time generated lookup table which maps the rank to char.
     *
     * A map generated at compile time where the key is the rank of the union
     * of all alternatives and the value is the corresponding char of that rank
     * and alternative.
     *
     */
    static constexpr std::array<char_type, value_size> value_to_char = []() constexpr
    {
        // Explicitly writing assign_rank_to_char within assign_value_to_char
        // causes this bug (g++-7 and g++-8):
        // https://gcc.gnu.org/bugzilla/show_bug.cgi?id=84684
        auto assign_rank_to_char = [](auto alternative, size_t rank) constexpr
        {
            return seqan3::to_char(seqan3::assign_rank(alternative, rank));
        };

        auto assign_value_to_char = [assign_rank_to_char] (auto alternative, auto & value_to_char, auto & value) constexpr
        {
            using alternative_t = std::decay_t<decltype(alternative)>;
            for (size_t i = 0u; i < alphabet_size_v<alternative_t>; ++i, ++value)
                value_to_char[value] = assign_rank_to_char(alternative, i);
        };

        unsigned value = 0u;
        std::array<char_type, value_size> value_to_char{};

        // initializer lists guarantee sequencing;
        // the following expression behaves as:
        // for(auto alternative: alternative_types)
        //    assign_value_to_char(alternative, value_to_char, value);
        ((assign_value_to_char(alternative_types{}, value_to_char, value)),...);

        return value_to_char;
    }();

    /*!\brief Compile-time generated lookup table which maps the char to rank.
     *
     * An map generated at compile time where the key is the char of one of the
     * alternatives and the value is the corresponding rank over all alternatives (by
     * conflict will default to the first).
     *
     */
    static constexpr std::array char_to_value = []() constexpr
    {
        constexpr size_t table_size = 1 << (sizeof(char_type) * 8);

        std::array<rank_type, table_size> char_to_value{};
        for (size_t i = 0u; i < value_to_char.size(); ++i)
        {
            using index_t = std::make_unsigned_t<char_type>;
            rank_type & old_entry = char_to_value[static_cast<index_t>(value_to_char[i])];
            bool is_new_entry = value_to_char[0] != value_to_char[i] && old_entry == 0;
            if (is_new_entry)
                old_entry = static_cast<rank_type>(i);
        }
        return char_to_value;
    }();

    //!\brief Converts an object of one of the given alternatives into the internal representation.
    //!\tparam index The position of `alternative_t` in the template pack `alternative_types`.
    //!\tparam alternative_t One of the alternative types.
    //!\param alternative The value of a alternative.
    template <size_t index, typename alternative_t>
    //!\cond
        requires holds_alternative<alternative_t>()
    //!\endcond
    static constexpr rank_type rank_by_index_(alternative_t const & alternative) noexcept
    {
        return partial_sum_sizes[index] + static_cast<rank_type>(seqan3::to_rank(alternative));
    }

    //!\brief Converts an object of one of the given alternatives into the internal representation.
    //!\details Finds the index of alternative_t in the given types.
    //!\tparam alternative_t One of the alternative types.
    //!\param alternative The value of a alternative.
    template <typename alternative_t>
    //!\cond
        requires holds_alternative<alternative_t>()
    //!\endcond
    static constexpr rank_type rank_by_type_(alternative_t const & alternative) noexcept
    {
        constexpr size_t index = meta::find_index<alternatives, alternative_t>::value;
        return rank_by_index_<index>(alternative);
    }
};

/*!\name Comparison operators
 * \relates union_composition
 * \brief Free function (in-)equality comparison operators that forward to member operators (for types != self).
 *\{
 */
template <typename lhs_t, typename ...alternative_types>
constexpr bool operator==(lhs_t const & lhs, union_composition<alternative_types...> const & rhs) noexcept
//!\cond
    requires detail::weakly_equality_comparable_by_members_with_concept<union_composition<alternative_types...>, lhs_t> &&
             !detail::weakly_equality_comparable_by_members_with_concept<lhs_t, union_composition<alternative_types...>>
//!\endcond
{
    return rhs == lhs;
}

template <typename lhs_t, typename ...alternative_types>
constexpr bool operator!=(lhs_t const & lhs, union_composition<alternative_types...> const & rhs) noexcept
//!\cond
    requires detail::weakly_equality_comparable_by_members_with_concept<union_composition<alternative_types...>, lhs_t> &&
             !detail::weakly_equality_comparable_by_members_with_concept<lhs_t, union_composition<alternative_types...>>
//!\endcond
{
    return rhs != lhs;
}
//!\}

} // namespace seqan3

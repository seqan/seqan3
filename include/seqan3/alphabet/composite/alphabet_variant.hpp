// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 * \author David Heller <david.heller AT fu-berlin.de>
 * \brief Provides seqan3::alphabet_variant.
 */

#pragma once

#include <algorithm>
#include <array>
#include <utility>
#include <cassert>
#include <variant>

#include <meta/meta.hpp>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/composite/detail.hpp>
#include <seqan3/alphabet/alphabet_base.hpp>
#include <seqan3/core/concept/core_language.hpp>
#include <seqan3/core/detail/int_types.hpp>
#include <seqan3/core/type_traits/pack.hpp>
#include <seqan3/core/type_traits/range.hpp>
#include <seqan3/core/type_traits/transformation_trait_or.hpp>
#include <seqan3/core/tuple_utility.hpp>
#include <seqan3/std/concepts>
#include <seqan3/std/type_traits>

namespace seqan3::detail
{

/*!\brief Evaluates to true if the one of the alternatives of the seqan3::alphabet_variant satisifes a compile-time
 *        predicate.
 * \tparam variant_t A specialisation of seqan3::alphabet_variant.
 * \tparam fun_t    A template template that takes target_t as argument and exposes an `invoke` member type that
 *                  evaluates some predicate and returns `std::true_type` or `std::false_type`.
 * \tparam target_t The type you wish query.
 * \ingroup composite
 *
 * \details
 *
 * To prevent recursive template and/or concept instantiation this call needs to be guarded against many exceptions.
 * See the source file for more details.
 */
// default is false
template <typename variant_t,
          template <typename> typename fun_t,
          typename target_t>
inline bool constexpr one_alternative_is = false;

//!\cond

// actual implementation
template <typename ... alternatives,
          template <typename> typename fun_t,
          typename target_t>
inline bool constexpr one_alternative_is<alphabet_variant<alternatives...>,
                                         fun_t,
                                         target_t>
 = !meta::empty<meta::find_if<meta::list<alternatives...>, fun_t<target_t>>>::value;

// guard against self
template <typename ... alternatives,
          template <typename> typename fun_t>
inline bool constexpr one_alternative_is<alphabet_variant<alternatives...>,
                                         fun_t,
                                         alphabet_variant<alternatives...>> = false;

// guard against types convertible to self without needing self's constructors
template <typename ... alternatives,
          template <typename> typename fun_t,
          typename target_t>
    requires ConvertibleToByMember<target_t, alphabet_variant<alternatives...>>
inline bool constexpr one_alternative_is<alphabet_variant<alternatives...>,
                                         fun_t,
                                         target_t> = false;

// guard against tuple composites that contain the variant somewhere (they can implicitly convert at source)
template <typename ... alternatives,
          template <typename> typename fun_t,
          typename target_t>
    requires AlphabetTupleBase<target_t> &&
             meta::in<detail::transformation_trait_or_t<recursive_tuple_components<target_t>, meta::list<>>,
                      alphabet_variant<alternatives...>>::value
inline bool constexpr one_alternative_is<alphabet_variant<alternatives...>,
                                         fun_t,
                                         target_t> = false;

// guard against alternatives
template <typename ... alternatives,
          template <typename> typename fun_t,
          typename target_t>
    requires type_in_pack_v<target_t, alternatives...>
inline bool constexpr one_alternative_is<alphabet_variant<alternatives...>,
                                         fun_t,
                                         target_t> = false;

// guard against alternatives (LHS and RHS switched)
template <typename ... alternatives,
          template <typename> typename fun_t,
          typename target_t>
    requires type_in_pack_v<target_t, alternatives...>
inline bool constexpr one_alternative_is<target_t,
                                         fun_t,
                                         alphabet_variant<alternatives...>> = false;

// guard against ranges and iterators over self to prevent recursive instantiation
template <typename ... alternatives,
          template <typename> typename fun_t,
          typename target_t>
    //NO, it's not possible to use the value_type type trait here
    requires requires { std::Same<typename target_t::value_type, alphabet_variant<alternatives...>>; }
inline bool constexpr one_alternative_is<alphabet_variant<alternatives...>,
                                         fun_t,
                                         target_t> = false;

// guard against pairs/tuples that *might* contain self to prevent recursive instantiation
// (applying TupleLike unfortunately does not work because it itself starts recursive instantiation)
template <typename ... alternatives,
          template <typename> typename fun_t,
          typename target_t>
    requires TupleSize<target_t> && !AlphabetTupleBase<target_t>
inline bool constexpr one_alternative_is<alphabet_variant<alternatives...>,
                                         fun_t,
                                         target_t> = false;

//!\endcond

} // namespace seqan3::detail

namespace seqan3
{

/*!\brief A combined alphabet that can hold values of either of its alternatives.
 * \ingroup composite
 * \if DEV
 * \tparam ...alternative_types Types of possible values (at least 2); all must model
 *                              seqan3::detail::WritableConstexprAlphabet, not be references and be unique.
 * \implements seqan3::detail::WritableConstexprAlphabet
 * \else
 * \tparam ...alternative_types Types of possible values (at least 2); all must model seqan3::WritableAlphabet,
 *                              must not be references and must be unique; all required functions for
 *                              seqan3::WritableAlphabet need to be callable in a `constexpr`-context.
 * \endif
 * \implements seqan3::WritableAlphabet
 * \implements seqan3::TriviallyCopyable
 * \implements seqan3::StandardLayout

 * \details
 *
 * The alphabet_variant represents the union of two or more alternative alphabets (e.g. the
 * four letter DNA alternative + the gap alternative). It behaves similar to a
 * [variant](https://en.cppreference.com/w/cpp/language/variant) or std::variant, but it preserves the
 * seqan3::Alphabet.
 *
 * Short description:
 *   * combines multiple different alphabets in an "either-or"-fashion;
 *   * is itself a seqan3::Alphabet;
 *   * its alphabet size is the sum of the individual sizes;
 *   * default initialises to the the first alternative's default (no empty state like std::variant);
 *   * constructible, assignable and (in-)equality-comparable with each alternative type and also all types that
 *     these are constructible/assignable/equality-comparable with;
 *   * only convertible to its alternatives through the member function convert_to() (which can throw!)
 *
 * ### Example
 *
 * \include test/snippet/alphabet/composite/alphabet_variant.cpp
 *
 * ### The `char` representation of an alphabet_variant
 *
 * Part of the seqan3::Alphabet concept requires that the alphabet_variant provides a char representation in addition
 * to the rank representation. For an object of seqan3::alphabet_variant, the `to_char()` member function will always
 * return the same character as if invoked on the respective alternative.
 * In contrast, the `assign_char()` member function might be ambiguous between the alternative alphabets in a variant.
 *
 * For example, assigning a '!' to seqan3::dna15 resolves to an object of rank 8 with char representation 'N' while
 * assigning '!' to seqan3::gap always resolves to rank 0, the gap symbol itself ('-'_gap).
 * We tackle this ambiguousness by **defaulting unknown characters to the representation of the first alternative**
 * (e.g. `alphabet_variant<dna15, gap>{}.assign_char('!')` resolves to rank 8, representing `N`_dna15).
 *
 * On the other hand, two alternative alphabets might have the same char representation (e.g if
 * you combine dna4 with dna5, 'A', 'C', 'G' and 'T' are ambiguous).
 * We tackle this ambiguousness by **always choosing the first valid char representation** (e.g.
 * `alphabet_variant<dna4, dna5>{}.assign_char('A')` resolves to rank 0, representing an `A`_dna4).
 *
 * To explicitly assign via the character representation of a specific alphabet,
 * assign to that type first and then assign to the variant, e.g.
 *
 * \include test/snippet/alphabet/composite/alphabet_variant_char_representation.cpp
 */
template <typename ...alternative_types>
//!\cond
    requires (detail::WritableConstexprAlphabet<alternative_types> && ...) &&
             (!std::is_reference_v<alternative_types> && ...) &&
             (sizeof...(alternative_types) >= 2)
             //TODO same char_type
//!\endcond
class alphabet_variant : public alphabet_base<alphabet_variant<alternative_types...>,
                                               (static_cast<size_t>(alphabet_size<alternative_types>) + ...),
                                               char> //TODO underlying char t

{
private:
    //!\brief The base type.
    using base_t = alphabet_base<alphabet_variant<alternative_types...>,
                                                   (static_cast<size_t>(alphabet_size<alternative_types>) + ...),
                                                   char>;
    //!\brief Befriend the base type.
    friend base_t;

    //!\brief A meta::list of the types of each alternative in the composite
    using alternatives = meta::list<alternative_types...>;

    static_assert(std::Same<alternatives, meta::unique<alternatives>>,
                  "All types in a alphabet_variant must be distinct.");

    using typename base_t::char_type;
    using typename base_t::rank_type;
public:
    using base_t::alphabet_size;
    using base_t::to_char;
    using base_t::to_rank;
    using base_t::assign_rank;

    /*!\brief Returns true if alternative_t is one of the given alternative types.
     * \tparam alternative_t The type to check.
     *
     * \include test/snippet/alphabet/composite/alphabet_variant_holds_alternative.cpp
     */
    template <typename alternative_t>
    static constexpr bool holds_alternative() noexcept
    {
        return detail::type_in_pack_v<alternative_t, alternative_types...>;
    }

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr alphabet_variant()                                     noexcept = default; //!< Defaulted.
    constexpr alphabet_variant(alphabet_variant const &)             noexcept = default; //!< Defaulted.
    constexpr alphabet_variant(alphabet_variant &&)                  noexcept = default; //!< Defaulted.
    constexpr alphabet_variant & operator=(alphabet_variant const &) noexcept = default; //!< Defaulted.
    constexpr alphabet_variant & operator=(alphabet_variant &&)      noexcept = default; //!< Defaulted.
    ~alphabet_variant()                                              noexcept = default; //!< Defaulted.

    /*!\brief Construction via the value of an alternative.
     * \tparam alternative_t One of the alternative types.
     * \param  alternative   The value of a alternative that should be assigned.
     *
     * \include test/snippet/alphabet/composite/alphabet_variant_value_construction.cpp
     */
    template <typename alternative_t>
    //!\cond
        requires holds_alternative<alternative_t>()
    //!\endcond
    constexpr alphabet_variant(alternative_t const & alternative) noexcept
    {
        assign_rank(rank_by_type_(alternative));
    }

    /*!\brief Construction via the value of a type that an alternative type is constructible from.
     * \tparam indirect_alternative_t A type that one of the alternative types is constructible from.
     * \param  rhs The value that should be assigned.
     *
     * \include test/snippet/alphabet/composite/alphabet_variant_conversion.cpp
     * \attention When selecting the alternative alphabet types which require only implicit conversion
     * or constructor calls, are preferred over those that require explicit ones.
     */
    template <typename indirect_alternative_t>
    //!\cond
        requires !detail::one_alternative_is<alphabet_variant,
                                             detail::implicitly_convertible_from,
                                             indirect_alternative_t> &&
                 detail::one_alternative_is<alphabet_variant,
                                            detail::constructible_from,
                                            indirect_alternative_t>
    //!\endcond
    constexpr alphabet_variant(indirect_alternative_t const & rhs) noexcept
    {
        assign_rank(rank_by_type_(meta::front<meta::find_if<alternatives,
                                                            detail::constructible_from<indirect_alternative_t>>>(rhs)));
    }

    //!\cond
    template <typename indirect_alternative_t>
        requires detail::one_alternative_is<alphabet_variant,
                                            detail::implicitly_convertible_from,
                                            indirect_alternative_t>
    constexpr alphabet_variant(indirect_alternative_t const & rhs) noexcept
    {
        assign_rank(
            rank_by_type_(
                meta::front<meta::find_if<alternatives,
                                          detail::implicitly_convertible_from<indirect_alternative_t>>>(rhs)));
    }
    //!\endcond

    /*!\brief Assignment via a value that one of the alternative types is assignable from.
     * \tparam indirect_alternative_t A type that one of the alternatives is assignable from.
     * \param  rhs The value of an alternative.
     *
     * \include test/snippet/alphabet/composite/alphabet_variant_subtype_construction.cpp
     */
    template <typename indirect_alternative_t>
    //!\cond
        requires !detail::one_alternative_is<alphabet_variant,
                                             detail::implicitly_convertible_from,
                                             indirect_alternative_t> &&             // constructor takes care
                 !detail::one_alternative_is<alphabet_variant,
                                             detail::constructible_from,
                                             indirect_alternative_t> &&             // constructor takes care
                 detail::one_alternative_is<alphabet_variant,
                                            detail::assignable_from,
                                            indirect_alternative_t>
    //!\endcond
    constexpr alphabet_variant & operator=(indirect_alternative_t const & rhs) noexcept
    {
        using alternative_t = meta::front<meta::find_if<alternatives, detail::assignable_from<indirect_alternative_t>>>;
        alternative_t alternative{};
        alternative = rhs;
        assign_rank(rank_by_type_(alternative));
        return *this;
    }
    //!\}

    /*!\name Conversion (by index)
     * \{
     */
    //!\brief Whether the variant alphabet currently holds a value of the given alternative.
    //!\tparam index Index of the alternative to check for.
    template <size_t index>
    constexpr bool is_alternative() const noexcept
    {
        static_assert(index < alphabet_size, "The alphabet_variant contains less alternatives than you are checking.");
        return (to_rank() >= partial_sum_sizes[index]) && (to_rank() < partial_sum_sizes[index + 1]);
    }

    /*!\brief Convert to the specified alphabet (throws if is_alternative() would be false).
     * \tparam index Index of the alternative to check for.
     * \throws std::bad_variant_access If the variant_alphabet currently holds the value of a different alternative.
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
     * \throws std::bad_variant_access If the variant_alphabet currently holds the value of a different alternative.
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

    /*!\name Comparison operators (against alternatives)
     * \brief Defines comparison against alternatives, e.g. `alphabet_variant<dna5, gap>{gap{}} == 'C'_dna5`. Only
     *        (in-)equality comparison is explicitly defined, because it would be difficult to argue about e.g.
     *        `alphabet_variant<dna5, gap>{gap{}} < 'C'_dna5`.
     * \{
     */
    template <typename alternative_t>
    constexpr bool operator==(alternative_t const rhs) const noexcept
        requires holds_alternative<alternative_t>()
    {
        return is_alternative<alternative_t>() && (convert_unsafely_to<alternative_t>() == rhs);
    }

    template <typename alternative_t>
    constexpr bool operator!=(alternative_t const rhs) const noexcept
        requires holds_alternative<alternative_t>()
    {
        return !operator==(rhs);
    }
    //!\}

    /*!\name Comparison operators (against indirect alternatives)
     * \brief Defines comparison against types that are comparable with alternatives, e.g.
     *        `alphabet_variant<dna5, gap>{'C'_dna5} == 'C'_rna5`. Only (in-)equality comparison is explicitly defined,
     *        because it would be difficult to argue about e.g.
     *        `alphabet_variant<dna5, gap>{gap{}} < 'C'_rna5`.
     * \{
     */
    template <typename indirect_alternative_type>
    constexpr bool operator==(indirect_alternative_type const rhs) const noexcept
    //!\cond
        requires detail::one_alternative_is<alphabet_variant,
                                            detail::weakly_equality_comparable_with,
                                            indirect_alternative_type>
    //!\endcond
    {
        using alternative_t =
            meta::front<meta::find_if<alternatives,
                                      detail::weakly_equality_comparable_with<indirect_alternative_type>>>;
        return is_alternative<alternative_t>() && (convert_unsafely_to<alternative_t>() == rhs);
    }

    template <typename indirect_alternative_type>
    constexpr bool operator!=(indirect_alternative_type const rhs) const noexcept
    //!\cond
        requires detail::one_alternative_is<alphabet_variant,
                                            detail::weakly_equality_comparable_with,
                                            indirect_alternative_type>
    //!\endcond
    {
        return !operator==(rhs);
    }
    //!\}

protected:
    //!\privatesection

    /*!\brief Implementation function for convert_to() and convert_unsafely_to().
     * \tparam index  Index of the alternative to convert to.
     * \tparam throws Whether to perform checks (and throw) or not.
     */
    template <size_t index, bool throws>
    constexpr auto convert_impl() const noexcept(!throws) -> meta::at_c<alternatives, index>
    {
        static_assert(index < alphabet_size, "The alphabet_variant contains less alternatives than you are checking.");
        using alternative_t = meta::at_c<alternatives, index>;

        if constexpr (throws)
        {
            if (!is_alternative<index>()) // [[unlikely]]
            {
                throw std::bad_variant_access{};
            }
        }

        using seqan3::assign_rank_to;
        return assign_rank_to(to_rank() - partial_sum_sizes[index], alternative_t{});
    }

    /*!\brief Compile-time generated lookup table which contains the partial
     * sum up to the position of each alternative.
     *
     * An array which contains the prefix sum over all
     * alternative_types::alphabet_size's.
     *
     */
    static constexpr std::array partial_sum_sizes = []() constexpr
    {
        constexpr size_t N = sizeof...(alternative_types) + 1;

        std::array<rank_type, N> partial_sum{0, seqan3::alphabet_size<alternative_types>...};
        for (size_t i = 1u; i < N; ++i)
            partial_sum[i] += partial_sum[i-1];

        return partial_sum;
    }();

    /*!\brief Compile-time generated lookup table which maps the rank to char.
     *
     * A map generated at compile time where the key is the rank of the variant
     * of all alternatives and the value is the corresponding char of that rank
     * and alternative.
     *
     */
    static constexpr std::array<char_type, alphabet_size> rank_to_char = []() constexpr
    {
        // Explicitly writing assign_rank_to_char within assign_rank_to_char
        // causes this bug (g++-7 and g++-8):
        // https://gcc.gnu.org/bugzilla/show_bug.cgi?id=84684
        auto assign_rank_to_char = [](auto alternative, size_t rank) constexpr
        {
            return seqan3::to_char(seqan3::assign_rank_to(rank, alternative));
        };

        auto assign_value_to_char = [assign_rank_to_char] (auto alternative, auto & value_to_char, auto & value) constexpr
        {
            using alternative_t = std::decay_t<decltype(alternative)>;
            for (size_t i = 0u; i < seqan3::alphabet_size<alternative_t>; ++i, ++value)
                value_to_char[value] = assign_rank_to_char(alternative, i);
        };

        unsigned value = 0u;
        std::array<char_type, alphabet_size> value_to_char{};

        // initializer lists guarantee sequencing;
        // the following expression behaves as:
        // for(auto alternative: alternative_types)
        //    assign_rank_to_char(alternative, rank_to_char, value);
        ((assign_value_to_char(alternative_types{}, value_to_char, value)),...);

        return value_to_char;
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

    /*!\brief Compile-time generated lookup table which maps the char to rank.
     *
     * An map generated at compile time where the key is the char of one of the
     * alternatives and the value is the corresponding rank over all alternatives (by
     * conflict will default to the first).
     *
     */
    static constexpr std::array char_to_rank = []() constexpr
    {
        constexpr size_t table_size = 1 << (sizeof(char_type) * 8);

        std::array<rank_type, table_size> char_to_rank{};

        for (size_t i = 0u; i < table_size; ++i)
        {
            char_type chr = static_cast<char_type>(i);
            bool there_was_no_valid_representation{true};

            meta::for_each(alternatives{}, [&] (auto && alt)
            {
                using alt_type = remove_cvref_t<decltype(alt)>;

                if (there_was_no_valid_representation && char_is_valid_for<alt_type>(chr))
                {
                    there_was_no_valid_representation = false;
                    char_to_rank[i] = rank_by_type_(assign_char_to(chr, alt_type{}));
                }
            });

            if (there_was_no_valid_representation)
                char_to_rank[i] = rank_by_type_(assign_char_to(chr, meta::front<alternatives>{}));
        }

        return char_to_rank;
    }();

    //!\brief Validate whether a character is valid in the by any of the combined alphabet.
    static constexpr bool char_is_valid(char_type const chr) noexcept
    {
        bool is_valid{false};

        meta::for_each(alternatives{}, [&] (auto && alt)
        {
            if (char_is_valid_for<remove_cvref_t<decltype(alt)>>(chr))
                is_valid = true;
        });

        return is_valid;
    }
};

/*!\name Comparison operators
 * \relates alphabet_variant
 * \brief Free function (in-)equality comparison operators that forward to member operators (for types != self).
 *\{
 */
template <typename lhs_t, typename ...alternative_types>
constexpr bool operator==(lhs_t const lhs, alphabet_variant<alternative_types...> const rhs) noexcept
//!\cond
    requires detail::WeaklyEqualityComparableByMembersWith<alphabet_variant<alternative_types...>, lhs_t> &&
             !detail::WeaklyEqualityComparableByMembersWith<lhs_t, alphabet_variant<alternative_types...>>
//!\endcond
{
    return rhs == lhs;
}

template <typename lhs_t, typename ...alternative_types>
constexpr bool operator!=(lhs_t const lhs, alphabet_variant<alternative_types...> const rhs) noexcept
//!\cond
    requires detail::WeaklyEqualityComparableByMembersWith<alphabet_variant<alternative_types...>, lhs_t> &&
             !detail::WeaklyEqualityComparableByMembersWith<lhs_t, alphabet_variant<alternative_types...>>
//!\endcond
{
    return rhs != lhs;
}
//!\}

} // namespace seqan3

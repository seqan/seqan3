// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \author David Heller <david.heller AT fu-berlin.de>
 * \brief Provides seqan3::alphabet_variant.
 */

#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <span>
#include <utility>
#include <variant>

#include <seqan3/alphabet/alphabet_base.hpp>
#include <seqan3/alphabet/composite/detail.hpp>
#include <seqan3/utility/type_list/traits.hpp>
#include <seqan3/utility/type_traits/lazy_conditional.hpp>

namespace seqan3::detail
{

//!\brief Helper to determine whether a type is a std::span.
template <typename t, size_t extent = std::dynamic_extent>
struct is_span : std::false_type
{};

//!\cond
template <typename t, size_t extent>
struct is_span<std::span<t, extent>> : std::true_type
{};
//!\endcond

//!\brief Prevents wrong instantiations of std::alphabet_variant's constructors.
template <typename other_t, typename... alternative_types>
inline constexpr bool variant_general_guard =
    (!std::same_as<other_t, alphabet_variant<alternative_types...>>)
    && (!std::is_base_of_v<alphabet_variant<alternative_types...>, other_t>)
    && (!(std::same_as<other_t, alternative_types> || ...))
    && (!list_traits::contains<alphabet_variant<alternative_types...>, recursive_required_types_t<other_t>>);

//!\brief Prevents wrong instantiations of std::alphabet_variant's comparison operators.
template <typename lhs_t, typename rhs_t, bool lhs_rhs_switched, typename... alternative_types>
inline constexpr bool variant_comparison_guard =
    (instantiate_if_v<lazy<weakly_equality_comparable_with_trait, rhs_t, alternative_types>,
                      (std::same_as<lhs_t, alphabet_variant<alternative_types...>>)
                          && (variant_general_guard<rhs_t, alternative_types...>)
                          && !(lhs_rhs_switched && is_type_specialisation_of_v<rhs_t, alphabet_variant>)>
     || ...);
} // namespace seqan3::detail

namespace seqan3
{

/*!\brief A combined alphabet that can hold values of either of its alternatives.
 * \ingroup alphabet_composite
 * \if DEV
 * \tparam ...alternative_types Types of possible values (at least 2); all must model
 *                              seqan3::detail::writable_constexpr_alphabet, std::regular and be unique.
 * \implements seqan3::detail::writable_constexpr_alphabet
 * \else
 * \tparam ...alternative_types Types of possible values (at least 2); all must model seqan3::writable_alphabet,
 *                              std::regular and must be unique; all required functions for
 *                              seqan3::writable_alphabet need to be callable in a `constexpr`-context.
 * \endif
 * \implements seqan3::writable_alphabet
 * \implements seqan3::trivially_copyable
 * \implements seqan3::standard_layout

 * \details
 *
 * The alphabet_variant represents the union of two or more alternative alphabets (e.g. the
 * four letter DNA alternative + the gap alternative). It behaves similar to a
 * [variant](https://en.cppreference.com/w/cpp/utility/variant) or std::variant, but it preserves the
 * seqan3::alphabet.
 *
 * Short description:
 *   * combines multiple different alphabets in an "either-or"-fashion;
 *   * is itself a seqan3::alphabet;
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
 * Part of the seqan3::alphabet concept requires that the alphabet_variant provides a char representation in addition
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
 *
 * \stableapi{Since version 3.1.}
 */
template <typename... alternative_types>
    requires (detail::writable_constexpr_alphabet<alternative_types> && ...) && (std::regular<alternative_types> && ...)
          && (sizeof...(alternative_types) >= 2)
class alphabet_variant :
    public alphabet_base<alphabet_variant<alternative_types...>,
                         (static_cast<size_t>(alphabet_size<alternative_types>) + ...),
                         char>
{
private:
    //!\brief The base type.
    using base_t = alphabet_base<alphabet_variant<alternative_types...>,
                                 (static_cast<size_t>(alphabet_size<alternative_types>) + ...),
                                 char>;

    static_assert((std::is_same_v<alphabet_char_t<alternative_types>, char> && ...),
                  "The alphabet_variant is currently only tested for alphabets with char_type char. "
                  "Contact us on GitHub if you have a different use case: https://github.com/seqan/seqan3 .");

    //!\brief Befriend the base type.
    friend base_t;

    //!\brief A seqan3::type_list of the types of each alternative in the composite
    using alternatives = seqan3::type_list<alternative_types...>;

    static_assert(((seqan3::list_traits::count<alternative_types, alternatives> == 1) && ... && true),
                  "All types in a alphabet_variant must be distinct.");

    using typename base_t::char_type;
    using typename base_t::rank_type;

public:
    using base_t::alphabet_size;
    using base_t::assign_rank;
    using base_t::to_char;
    using base_t::to_rank;

    /*!\brief Expose the alternative types to concept checks in metaprogramming.
    * \private
    * \details
    * \noapi
    */
    using seqan3_required_types = type_list<alternative_types...>;
    /*!\brief Expose the recursive alternative types to concept checks in metaprogramming.
    * \private
    * \details
    * \noapi
    */
    using seqan3_recursive_required_types = list_traits::concat<
        seqan3_required_types,
        detail::transformation_trait_or_t<detail::recursive_required_types<alternative_types>, type_list<>>...>;

    /*!\brief Returns true if alternative_t is one of the given alternative types.
     * \tparam alternative_t The type to check.
     *
     * \include test/snippet/alphabet/composite/alphabet_variant_is_alternative.cpp
     *
     * \stableapi{Since version 3.1.}
     */
    template <typename alternative_t>
    static constexpr bool is_alternative() noexcept
    {
        return seqan3::pack_traits::contains<alternative_t, alternative_types...>;
    }

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr alphabet_variant() noexcept = default;                                     //!< Defaulted.
    constexpr alphabet_variant(alphabet_variant const &) noexcept = default;             //!< Defaulted.
    constexpr alphabet_variant(alphabet_variant &&) noexcept = default;                  //!< Defaulted.
    constexpr alphabet_variant & operator=(alphabet_variant const &) noexcept = default; //!< Defaulted.
    constexpr alphabet_variant & operator=(alphabet_variant &&) noexcept = default;      //!< Defaulted.
    ~alphabet_variant() noexcept = default;                                              //!< Defaulted.

    /*!\brief Construction via the value of an alternative.
     * \tparam alternative_t One of the alternative types.
     * \param  alternative   The value of a alternative that should be assigned.
     *
     * \include test/snippet/alphabet/composite/alphabet_variant_value_construction.cpp
     *
     * \stableapi{Since version 3.1.}
     */
    template <typename alternative_t>
        requires (!std::same_as<alternative_t, alphabet_variant>) && (!detail::is_span<alternative_t>::value)
              && (!std::is_base_of_v<alphabet_variant, alternative_t>)
              && (!list_traits::contains<
                  alphabet_variant,
                  detail::transformation_trait_or_t<detail::recursive_required_types<alternative_t>, type_list<>>>)
              && (is_alternative<alternative_t>())
    constexpr alphabet_variant(alternative_t const alternative) noexcept
    {
        assign_rank(rank_by_type_(alternative));
    }

    /*!\brief Constructor for arguments implicitly convertible to an alternative.
     * \tparam indirect_alternative_t A type that is implicitly convertible to an alternative type.
     * \param  rhs The value that should be assigned.
     *
     * \details
     *
     * This constructor is preferred over the explicit version.
     *
     * ### Example
     *
     * \include test/snippet/alphabet/composite/alphabet_variant_conversion.cpp
     *
     *   * seqan3::dna4 and seqan3::rna4 are implicitly convertible to each other so the variant accepts either.
     *   * Construction via `{}` considers implicit and explicit conversions.
     *   * Construction via `=` considers only implicit conversions (but that is sufficient here).
     *
     * \experimentalapi{Experimental since version 3.1.}
     */
    template <typename indirect_alternative_t>
        requires (!detail::is_span<indirect_alternative_t>::value)
              && ((
                  detail::instantiate_if_v<detail::lazy<std::is_convertible, indirect_alternative_t, alternative_types>,
                                           detail::variant_general_guard<indirect_alternative_t, alternative_types...>>
                  || ...))
    constexpr alphabet_variant(indirect_alternative_t const rhs) noexcept
    {
        using alternative_predicate = detail::implicitly_convertible_from<indirect_alternative_t>;
        constexpr auto alternative_position =
            seqan3::list_traits::find_if<alternative_predicate::template invoke, alternatives>;
        using alternative_t = seqan3::list_traits::at<alternative_position, alternatives>;
        assign_rank(rank_by_type_(alternative_t(rhs)));
    }

    /*!\brief Constructor for arguments explicitly (but not implicitly) convertible to an alternative.
     * \tparam indirect_alternative_t A type that is explicitly (but not implicitly) convertible to an alternative type.
     * \param  rhs The value that should be assigned.
     *
     * \details
     *
     * ### Example
     *
     * \include test/snippet/alphabet/composite/alphabet_variant_conversion_explicit.cpp
     *
     *   * seqan3::dna4 and seqan3::dna5 are not implicitly convertible to each other, only explicitly.
     *   * Construction via `{}` considers implicit and explicit conversions so this works.
     *   * Construction via `=` considers only implicit conversions so it does not work.
     *
     * \experimentalapi{Experimental since version 3.1.}
     */
    template <typename indirect_alternative_t>
        requires (
            (!(detail::instantiate_if_v<detail::lazy<std::is_convertible, indirect_alternative_t, alternative_types>,
                                        detail::variant_general_guard<indirect_alternative_t, alternative_types...>>
               || ...))
            && (detail::instantiate_if_v<detail::lazy<std::is_constructible, alternative_types, indirect_alternative_t>,
                                         detail::variant_general_guard<indirect_alternative_t, alternative_types...>>
                || ...))
    constexpr explicit alphabet_variant(indirect_alternative_t const rhs) noexcept
    {
        using alternative_predicate = detail::constructible_from<indirect_alternative_t>;
        constexpr auto alternative_position =
            seqan3::list_traits::find_if<alternative_predicate::template invoke, alternatives>;
        using alternative_t = seqan3::list_traits::at<alternative_position, alternatives>;
        assign_rank(rank_by_type_(alternative_t(rhs)));
    }

    /*!\brief Assignment for arguments assignable to an alternative.
     * \tparam indirect_alternative_t A type that one of the alternatives is assignable from.
     * \param  rhs The value of an alternative.
     *
     * \details
     *
     * Most assignments happen through implicit conversion and the default assignment operator. This assignment operator
     * is for the rest.
     *
     * \experimentalapi{Experimental since version 3.1.}
     */
    template <typename indirect_alternative_t>
        requires (detail::variant_general_guard<indirect_alternative_t, alternative_types...>
                  && (weakly_assignable_from<alternative_types, indirect_alternative_t> || ...))
    constexpr alphabet_variant & operator=(indirect_alternative_t const & rhs) noexcept
    {
        using alternative_predicate = detail::assignable_from<indirect_alternative_t>;
        constexpr auto alternative_position =
            seqan3::list_traits::find_if<alternative_predicate::template invoke, alternatives>;
        using alternative_t = seqan3::list_traits::at<alternative_position, alternatives>;
        alternative_t alternative{};
        alternative = rhs;
        assign_rank(rank_by_type_(alternative));
        return *this;
    }
    //!\}

    /*!\name Conversion (by index)
     * \{
     */
    /*!\brief Whether the variant alphabet currently holds a value of the given alternative.
     * \tparam index Index of the alternative to check for.
     * \sa std::variant::holds_alternative
     * \details
     * \experimentalapi{Experimental since version 3.1.}
     */
    template <size_t index>
    constexpr bool holds_alternative() const noexcept
    {
        static_assert(index < alphabet_size, "The alphabet_variant contains less alternatives than you are checking.");
        return (to_rank() >= partial_sum_sizes[index]) && (to_rank() < partial_sum_sizes[index + 1]);
    }

    /*!\brief Convert to the specified alphabet (throws if holds_alternative() would be false).
     * \tparam index Index of the alternative to check for.
     * \throws std::bad_variant_access If the variant_alphabet currently holds the value of a different alternative.
     * \details
     * \experimentalapi{Experimental since version 3.1.}
     */
    template <size_t index>
    constexpr auto convert_to() const
    {
        return convert_impl<index, true>();
    }

    /*!\brief Convert to the specified alphabet (**undefined behaviour** if holds_alternative() would be false).
     * \tparam index Index of the alternative to check for.
     * \details
     * \experimentalapi{Experimental since version 3.1.}
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
    /*!\copybrief holds_alternative()
     * \tparam alternative_t The type of the alternative that you wish to check for.
     * \details
     * \stableapi{Since version 3.1.}
     */
    template <typename alternative_t>
    constexpr bool holds_alternative() const noexcept
        requires (is_alternative<alternative_t>())
    {
        constexpr size_t index = seqan3::list_traits::find<alternative_t, alternatives>;
        return holds_alternative<index>();
    }

    /*!\copybrief convert_to()
     * \tparam alternative_t The type of the alternative that you wish to check for.
     * \throws std::bad_variant_access If the variant_alphabet currently holds the value of a different alternative.
     * \details
     * \experimentalapi{Experimental since version 3.1.}
     */
    template <typename alternative_t>
    constexpr alternative_t convert_to() const
        requires (is_alternative<alternative_t>())
    {
        constexpr size_t index = seqan3::list_traits::find<alternative_t, alternatives>;
        return convert_impl<index, true>();
    }

    /*!\copybrief convert_unsafely_to()
     * \tparam alternative_t The type of the alternative that you wish to check for.
     * \details
     * \experimentalapi{Experimental since version 3.1.}
     */
    template <typename alternative_t>
    constexpr alternative_t convert_unsafely_to() const noexcept
        requires (is_alternative<alternative_t>())
    {
        constexpr size_t index = seqan3::list_traits::find<alternative_t, alternatives>;
        return convert_impl<index, false>();
    }
    //!\}

    /*!\name Comparison operators (against indirect alternatives)
     * \brief Defines comparison against types that are not subject to implicit construction/conversion but are
     *        comparable against alternatives, e.g. `alphabet_variant<seqan3::rna4, seqan3::gap>` vs
     *        `alphabet_variant<seqan3::dna4, seqan3::gap>`. Only (in-)equality comparison is defined as reasoning
     *        about order of variants is inherently difficult.
     * \{
     */
    /*!\brief (In-)Equality comparison against types comparable with alternatives but not convertible to the variant.
     * \tparam alphabet_variant_t The type of the variant; given as template parameter to prevent conversion.
     * \tparam indirect_alternative_type Must be comparable with an alternative's type.
     * \param lhs Left-hand-side of comparison.
     * \param rhs Right-hand-side of comparison.
     * \returns `true` or `false`.
     *
     * \details
     *
     * To determine (in-)equality, it is first deduced which alternative the argument is comparable with.
     * It is then checked if the variant currently is in that alternative's state and if yes whether the values compare
     * to `true`; else `false` is returned.
     * \details
     * \stableapi{Since version 3.1.}
     */
    template <typename alphabet_variant_t, typename indirect_alternative_type>
        requires (!detail::is_span<alphabet_variant_t>::value)
              && (detail::variant_comparison_guard<alphabet_variant_t,
                                                   indirect_alternative_type,
                                                   false,
                                                   alternative_types...>)
    friend constexpr bool operator==(alphabet_variant_t const lhs, indirect_alternative_type const rhs) noexcept
    {
        using alternative_predicate = detail::weakly_equality_comparable_with_<indirect_alternative_type>;
        constexpr auto alternative_position =
            seqan3::list_traits::find_if<alternative_predicate::template invoke, alternatives>;
        using alternative_t = seqan3::list_traits::at<alternative_position, alternatives>;
        return lhs.template holds_alternative<alternative_t>()
            && (lhs.template convert_unsafely_to<alternative_t>() == rhs);
    }

    //!\copydoc operator==(alphabet_variant_t const lhs, indirect_alternative_type const rhs)
    template <typename alphabet_variant_t, typename indirect_alternative_type>
        requires (!detail::is_span<alphabet_variant_t>::value)
              && (detail::variant_comparison_guard<alphabet_variant_t,
                                                   indirect_alternative_type,
                                                   false,
                                                   alternative_types...>)
    friend constexpr bool operator!=(alphabet_variant_t const lhs, indirect_alternative_type const rhs) noexcept
    {
        return !(lhs == rhs);
    }

    //!\copydoc operator==(alphabet_variant_t const lhs, indirect_alternative_type const rhs)
    template <typename alphabet_variant_t, typename indirect_alternative_type, typename = void>
        requires (!detail::is_span<alphabet_variant_t>::value)
              && (detail::variant_comparison_guard<alphabet_variant_t,
                                                   indirect_alternative_type,
                                                   true,
                                                   alternative_types...>)
    friend constexpr bool operator==(indirect_alternative_type const lhs, alphabet_variant_t const rhs) noexcept
    {
        return rhs == lhs;
    }

    //!\copydoc operator==(alphabet_variant_t const lhs, indirect_alternative_type const rhs)
    template <typename alphabet_variant_t, typename indirect_alternative_type, typename = void>
        requires (!detail::is_span<alphabet_variant_t>::value)
              && (detail::variant_comparison_guard<alphabet_variant_t,
                                                   indirect_alternative_type,
                                                   true,
                                                   alternative_types...>)
    friend constexpr bool operator!=(indirect_alternative_type const lhs, alphabet_variant_t const rhs) noexcept
    {
        return rhs != lhs;
    }
    //!\}

    /*!\brief Validate whether a character is valid in the combined alphabet.
     *
     * \experimentalapi{Experimental since version 3.1.}
     */
    static constexpr bool char_is_valid(char_type const chr) noexcept
    {
        using index_t = std::make_unsigned_t<char_type>;
        return first_valid_char_table[static_cast<index_t>(chr)] < sizeof...(alternative_types);
    }

protected:
    //!\privatesection

    /*!\brief Implementation function for convert_to() and convert_unsafely_to().
     * \tparam index  Index of the alternative to convert to.
     * \tparam throws Whether to perform checks (and throw) or not.
     */
    template <size_t index, bool throws>
    constexpr auto convert_impl() const noexcept(!throws) -> seqan3::list_traits::at<index, alternatives>
    {
        static_assert(index < alphabet_size, "The alphabet_variant contains less alternatives than you are checking.");
        using alternative_t = seqan3::list_traits::at<index, alternatives>;

        if constexpr (throws)
        {
            if (!holds_alternative<index>()) // [[unlikely]]
            {
                throw std::bad_variant_access{};
            }
        }

        return seqan3::assign_rank_to(to_rank() - partial_sum_sizes[index], alternative_t{});
    }

    /*!\brief Converts an object of one of the given alternatives into the internal representation.
     * \tparam index The position of `alternative_t` in the template pack `alternative_types`.
     * \tparam alternative_t One of the alternative types.
     * \param alternative The value of a alternative.
     */
    template <size_t index, typename alternative_t>
        requires (is_alternative<alternative_t>())
    static constexpr rank_type rank_by_index_(alternative_t const & alternative) noexcept
    {
        return partial_sum_sizes[index] + static_cast<rank_type>(seqan3::to_rank(alternative));
    }

    /*!\brief Converts an object of one of the given alternatives into the internal representation.
     * \details Finds the index of alternative_t in the given types.
     * \tparam alternative_t One of the alternative types.
     * \param alternative The value of a alternative.
     */
    template <typename alternative_t>
        requires (is_alternative<alternative_t>())
    static constexpr rank_type rank_by_type_(alternative_t const & alternative) noexcept
    {
        constexpr size_t index = seqan3::list_traits::find<alternative_t, alternatives>;
        return rank_by_index_<index>(alternative);
    }

    /*!\brief Compile-time generated lookup table which maps the rank to char.
     *
     * A map generated at compile time where the key is the rank of the variant
     * of all alternatives and the value is the corresponding char of that rank
     * and alternative.
     */
    static constexpr char_type rank_to_char(rank_type const rank)
    {
        return rank_to_char_table[rank];
    }

    /*!\brief Compile-time generated lookup table which maps the char to rank.
        *
        * A map generated at compile time where the key is the char of one of the
        * alternatives and the value is the corresponding rank over all alternatives (by
        * conflict will default to the first).
        */
    static constexpr rank_type char_to_rank(char_type const chr)
    {
        using index_t = std::make_unsigned_t<char_type>;
        return char_to_rank_table[static_cast<index_t>(chr)];
    }

    /*!\brief Compile-time generated lookup table which contains the partial
     * sum up to the position of each alternative.
     *
     * An array which contains the prefix sum over all
     * alternative_types::alphabet_size's.
     *
     */
    static constexpr std::array<rank_type, sizeof...(alternative_types) + 1> partial_sum_sizes{
        []() constexpr
        {
            constexpr size_t N = sizeof...(alternative_types) + 1;

            std::array<rank_type, N> partial_sum{0, seqan3::alphabet_size<alternative_types>...};

            for (size_t i = 1u; i < N; ++i)
                partial_sum[i] += partial_sum[i - 1];

            return partial_sum;
        }()};

    //!\copydoc seqan3::alphabet_variant::rank_to_char
    static constexpr std::array<char_type, alphabet_size> rank_to_char_table{
        []() constexpr
        {
            auto assign_value_to_char = [](auto alternative, auto & value_to_char, auto & value) constexpr
            {
                using alternative_t = std::decay_t<decltype(alternative)>;
                for (size_t i = 0u; i < seqan3::alphabet_size<alternative_t>; ++i, ++value)
                    value_to_char[value] = seqan3::to_char(seqan3::assign_rank_to(i, alternative));
            };

            size_t value{};
            std::array<char_type, alphabet_size> value_to_char{};

            // initializer lists guarantee sequencing;
            // the following expression behaves as:
            // for(auto alternative: alternative_types)
            //    assign_rank_to_char(alternative, rank_to_char, value);
            ((assign_value_to_char(alternative_types{}, value_to_char, value)), ...);

            return value_to_char;
        }()};

    /*!\brief Compile-time generated lookup table which maps the char to the index of the first alphabet that fulfils
     *        char_is_valid_for.
     */
    static constexpr auto first_valid_char_table{
        []() constexpr
        {
            constexpr size_t alternative_size = sizeof...(alternative_types);
            constexpr size_t table_size = detail::size_in_values_v<char_type>;
            using first_alphabet_t = detail::min_viable_uint_t<alternative_size>;

            std::array<first_alphabet_t, table_size> lookup_table{};

            for (size_t i = 0u; i < table_size; ++i)
            {
                char_type chr = static_cast<char_type>(i);

                std::array<bool, alternative_size> valid_chars{char_is_valid_for<alternative_types>(chr)...};

                auto found_it = std::find(valid_chars.begin(), valid_chars.end(), true);
                lookup_table[i] = found_it - valid_chars.begin();
            }

            return lookup_table;
        }()};

    //!\copydoc seqan3::alphabet_variant::char_to_rank
    static constexpr std::array<rank_type, detail::size_in_values_v<char_type>> char_to_rank_table{
        []() constexpr
        {
            constexpr size_t alternative_size = sizeof...(alternative_types);
            constexpr size_t table_size = detail::size_in_values_v<char_type>;

            std::array<rank_type, table_size> char_to_rank{};

            for (size_t i = 0u; i < table_size; ++i)
            {
                char_type chr = static_cast<char_type>(i);

                std::array<rank_type, alternative_size> ranks{
                    rank_by_type_(assign_char_to(chr, alternative_types{}))...};

                // if no char_is_valid_for any alternative, use the rank of the first alternative
                char_to_rank[i] = first_valid_char_table[i] < alternative_size ? ranks[first_valid_char_table[i]] : 0;
            }

            return char_to_rank;
        }()};
};

} // namespace seqan3

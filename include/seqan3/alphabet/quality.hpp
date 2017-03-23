#pragma once

#include <iostream>
#include <string>

#include "alphabet.hpp"

namespace seqan3
{

namespace detail
{

    template <typename q>
    concept bool internal_quality_concept = requires (q quality)
    {
        requires internal_alphabet_concept<q>;
        typename q::phred_type;


    { quality.to_phred() } -> typename q::phred_type;

    { quality.from_phred(0) } -> q;

    };
} // namespace seqan3::detail

template <typename alphabet_type>
requires detail::internal_quality_concept<alphabet_type>
struct underlying_phred
{
    using type = typename alphabet_type::phred_type;
};

template <typename alphabet_type>
using underlying_phred_t = typename underlying_phred<alphabet_type>::type;

template <typename alphabet_type>
requires detail::internal_quality_concept<alphabet_type>
constexpr alphabet_type from_phred(alphabet_type & c, char const in)
{
    return c.from_phred(in);
}

template <typename alphabet_type>
requires detail::internal_quality_concept<alphabet_type>
constexpr alphabet_char_t<alphabet_type> to_phred(alphabet_type const & c)
{
    return c.to_phred();
}

// ------------------------------------------------------------------
// concept
// ------------------------------------------------------------------

template<typename q>
concept bool quality_concept = requires(q quality)
{
    // requires fulfillment of alphabet concept
    requires alphabet_concept<q>;

    // conversion functions from t -> c and vice versa
    { from_phred(quality, typename q::integral_type{}) } -> q;
    { to_phred(quality) } -> const typename q::phred_type;
    typename underlying_phred<q>::type;
};

}
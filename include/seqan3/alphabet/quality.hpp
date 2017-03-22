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
constexpr alphabet_type from_char(alphabet_type & c, char const in)
{
    return c.from_phred(in);
}

template <typename alphabet_type>
requires detail::internal_quality_concept<alphabet_type>
constexpr alphabet_char_t<alphabet_type> to_char(alphabet_type const & c)
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

    // static data members of type value_type and char_type
    requires std::is_same_v<typename q::integral_type, decltype(q::value)>;

    // offsets for char and integer encoding
    q::offset_phred;
    q::offset_char;
    requires std::is_same_v<typename q::phred_type, std::decay_t<decltype(q::offset_phred)>>;
    requires std::is_same_v<typename q::char_type, std::decay_t<decltype(q::offset_char)>>;

    // conversion functions from t -> c and vice versa
    { quality.from_phred(typename q::integral_type{}) } -> q;
    { quality.to_phred() } -> const typename q::phred_type;


};

}
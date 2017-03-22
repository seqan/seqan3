#pragma once

#include <iostream>
#include <string>

#include "alphabet.hpp"

namespace seqan3 {

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
    requires std::is_same_v<typename q::phred_type, decltype(q::offset_phred)>;
    requires std::is_same_v<typename q::char_type, decltype(q::offset_char)>;

    // conversion functions from t -> c and vice versa
    { quality.to_phred() } -> typename q::phred_type;
    { quality.from_phred(typename q::integral_type{}) } -> q;


};

}
#pragma once

#include <optional>

#include "../alphabet.hpp"

namespace seqan3
{

namespace detail
{
// auto constexpr uint_min(uint64_t value)
// {
//     if constexpr(value <= std::numeric_limits<uint8_t>::max())
//         return uint8_t{};
//     else
//         return uint16_t{};
// }

}
template <typename underlying_t>
    requires alphabet_concept<underlying_t>
struct gaps
{
    /* types */
    using c_type = bool;

    /* member */
    c_type value;

    /* static member */
    static char gap_symbol{'-'};

    static constexpr c_type value_size{underlying_type::value_size + 1};

    /* public member functions */
    constexpr char to_char() const
    {
        return value_to_char[value];
    }

    constexpr c_type to_integral() const
    {
        return value;
    }

    constexpr gapped_alphabet from_char(char const in)
    {
        value = char_to_value[in];
        return *this;
    }

    constexpr gapped_alphabet from_integral(c_type const in)
    {
        value = in;
        return *this;
    }

    // conversion tables
    static constexpr uint8_t value_size{2};

    static constexpr char value_to_char[value_size]
    {
        ' ',
        '-'
    };

    static constexpr c_type char_to_value[256]
    {
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, //0
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, //1
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 1, 0, 0, 0, //2
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, //3
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, //4
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, //5
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, //6
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, //7
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, //8
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, //9
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, //10
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, //11
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, //12
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, //13
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, //14
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0  //15
    };
};

// static_assert(alphabet_concept<gapped_alphabet<dna4>>);

}

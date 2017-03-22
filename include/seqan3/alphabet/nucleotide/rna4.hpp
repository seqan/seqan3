#pragma once

#include "../alphabet.hpp"

namespace seqan3
{


struct rna4 : dna4
{

    // strictly typed enum, unfortunately with scope
    enum struct c_type : uint8_t

    // the value
    c_type value;

    constexpr rna4 & operator =(c_type const c)
    {
        value = c;
    }
    constexpr rna4 from_char(char const c)
    {
        value = char_to_value[c];
        return *this;
    }
    constexpr rna4 from_integral(std::underlying_type_t<c_type> const c)
    {
        value = static_cast<c_type>(c);
        return *this;
    }
    static constexpr uint8_t value_size{4};

    static constexpr char value_to_char[value_size]
    {
        'A',
        'C',
        'G',
        'U'
    };

   
};

// shall fulfill Alphabet concept
static_assert(alphabet_concept<rna4>);
static_assert(rna4{rna4::U} == rna4{});
static_assert(rna4{rna4::U} == rna4::U);
// static_assert(rna4{'U'} == 'U');
static_assert(static_cast<char>(rna4{rna4::C}) == 'C');
static_assert(rna4{rna4::A} < rna4{rna4::C});

}

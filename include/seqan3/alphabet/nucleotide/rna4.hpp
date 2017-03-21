#pragma once

#include "../alphabet.hpp"

namespace seqan3
{


struct rna4 : dna4
{
    // strictly typed enum, unfortunately with scope
    enum struct c_type : uint8_t
    {
        A,
        C,
        G,
        T,
        U = T,
        UNKNOWN = A
    };

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
        'T',
        'U'
    };

    static constexpr c_type char_to_value[256]
    {
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN,
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, //0
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN,
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, //1
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN,
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, //2
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN,
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, //3
        //          , A,            B,            C,            D,            E,            F,            G,
        c_type::UNKNOWN, c_type::A,       c_type::UNKNOWN, c_type::C,       c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::G,
        // H,         I,            J,            K,            L,            M,            N,            O,
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, //4
        // P,         Q,            R,            S,            T,            U,            V,            W,
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::T,       c_type::T,       c_type::UNKNOWN, c_type::UNKNOWN,
        // X,         Y,            Z,
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, //5
        //          , a,            b,            c,            d,            e,            f,            g,
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::C,       c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::G,
        // h,         i,            j,            k,            l,            m,            n,            o,
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, //6
        // P,         Q,            R,            S,            T,            U,            V,            W,
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::T,       c_type::T,       c_type::UNKNOWN, c_type::UNKNOWN,
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, //7
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN,
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, //8
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN,
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, //9
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN,
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, //10
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN,
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, //11
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN,
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, //12
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN,
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, //13
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN,
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, //14
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN,
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN  //15
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

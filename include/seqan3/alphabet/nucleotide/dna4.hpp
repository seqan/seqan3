#pragma once

#include "../alphabet.hpp"

namespace seqan3
{

struct dna4
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

    // import into local scope
    static constexpr c_type A{c_type::A};
    static constexpr c_type C{c_type::C};
    static constexpr c_type G{c_type::G};
    static constexpr c_type T{c_type::T};
    static constexpr c_type U{c_type::U};
    static constexpr c_type UNKNOWN{c_type::UNKNOWN};

    // the value
    c_type value;

    // implicit compatibility to inner_type
    constexpr dna4 & operator =(c_type const c)
    {
        value = c;
    }
    constexpr operator c_type() const
    {
        return value;
    }

    // explicit compatibility to char
    explicit constexpr operator char() const
    {
        return to_char();
    }
    constexpr char to_char() const
    {
        return value_to_char[static_cast<uint8_t>(value)];
    }

    constexpr dna4 from_char(char const c)
    {
        value = char_to_value[c];
        return *this;
    }

    // explicit compatibility to integral
    constexpr char to_integral() const
    {
        return static_cast<std::underlying_type_t<c_type>>(value);
    }

    constexpr dna4 from_integral(std::underlying_type_t<c_type> const c)
    {
        value = static_cast<c_type>(c);
        return *this;
    }

    // conversion tables
    static constexpr uint8_t value_size{4};

    static constexpr char value_to_char[value_size]
    {
        'A',
        'C',
        'G',
        'T'
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
static_assert(alphabet_concept<dna4>);
static_assert(dna4{dna4::A} == dna4{});
static_assert(dna4{dna4::A} == dna4::A);
// static_assert(dna4{'A'} == 'A');
static_assert(static_cast<char>(dna4{dna4::C}) == 'C');
static_assert(dna4{dna4::A} < dna4{dna4::C});

}

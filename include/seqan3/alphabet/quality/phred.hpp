#pragma once

#include "../alphabet.hpp"

#include <tuple>

// assume sequential mapping w.r.t. ascii alphabet, only (char_start, val_start) and (char_end, val_end) differ
namespace seqan3
{
    // 0 -> '!' ... 41 -> 'J'
    struct illumina18
    {
        // strictly typed enum, unfortunately with scope
        using c_type = uint8_t;
        char offset = '!';
        uint8_t size = 42;

        // the value
        c_type value;

        // implicit compatibility to inner_type
        constexpr illumina18 & operator =(c_type const c)
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

        //
        constexpr char to_char() const
        {
            return this->value + this->offset;
        }

        constexpr illumina18 from_char(char const c)
        {
            value = c - '!';
            return *this;
        }

        // explicit compatibility to integral
        constexpr char to_integral() const
        {
            return value;
        }

        constexpr illumina18 from_integral(uint8_t const c)
        {
            value = c;
            return *this;
        }

        // conversion tables
        static constexpr uint8_t value_size{};

    };

}

#pragma once

#include "gaps.hpp"

namespace seqan3
{

template <typename underlying_t>
    requires alphabet_concept<underlying_t>
struct gapped_alphabet : public compound_alphabet<underlying_t, gaps>
{
    /* public member functions */
    constexpr char to_char() const
    {
        return std::get<1>(*this).to_integral() ? std::get<1>(*this).to_char() : std::get<0>(*this).to_char();
    }

    constexpr char to_integral() const
    {
        return std::get<1>(*this).to_integral() ? std::get<1>(*this).to_integral() : std::get<0>(*this).to_integral();
    }

    constexpr bool is_gap() const
    {
        return value == value_size - 1;
    }

    constexpr gapped_alphabet from_char(char const in)
    {
        if (in == gap_symbol)
            value = value_size - 1;
        else
            value = underlying_type{}.from_char(in).to_integral();
        return *this;
    }

    constexpr gapped_alphabet from_integral(c_type const in)
    {
        if (in == value_size - 1)
            value.reset();
        else
            value = underlying_type{}.from_integral(in);
        return *this;
    }

    constexpr gapped_alphabet set_gap()
    {
        value.reset();
        return *this;
    }
};

// static_assert(alphabet_concept<gapped_alphabet<dna4>>);

}

#pragma once

template <typename ...alph_ts>
    requires (alphabet_concept<alph_ts> && ...)
struct compound_alphabet : public std::tuple<alph_ts...>
{
    constexpr char to_char() const
    {
        return std::get<0>(*this).to_char();
    }

    constexpr auto to_integral() const
    {
        return std::get<0>(*this).to_integral();
    }

    constexpr auto from_char(char const c)
    {
        return std::get<0>(*this).from_char(c);
    }

    //TODO auto select integral type
    constexpr auto from_integral(uint8_t const i)
    {
        // TODO return proper value
        return (std::get<0>(*this).from_integral(i) | ;
    }

    // todo autodetect type
    static constexpr uint8_t value_size = (alph_ts::value_size * ...);
};

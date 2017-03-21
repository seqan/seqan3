#pragma once

#include <iostream>
#include <string>
#include <vector>

#include "../alphabet.hpp"
#include "../alphabet_container.hpp"
#include "dna4.hpp"

// ------------------------------------------------------------------
// containers
// -----------------------------------------------------------------

namespace seqan3
{

using dna4_vector = std::vector<dna4>;

using dna4_string = std::basic_string<dna4, std::char_traits<dna4>>;

} // namespace seqan3

// ------------------------------------------------------------------
// literals
// -----------------------------------------------------------------

namespace seqan3::literal
{

inline dna4_vector operator "" _dna4(const char * s, std::size_t n)
{
    dna4_vector r;
    r.resize(n);

    std::transform(s, s + n, r.begin(), [] (const char & c)
    {
        return dna4{dna4::char_to_value[c]};
    });

    return r;
}

inline dna4_string operator "" _dna4s(const char * s, std::size_t n)
{
    dna4_string r;
    r.resize(n);

    std::transform(s, s + n, r.begin(), [] (const char & c)
    {
        return dna4{dna4::char_to_value[c]};
    });

    return r;
}

} // namespace seqan3::literal


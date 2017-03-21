#pragma once

#include <iostream>
#include <string>
#include <vector>

#include "../alphabet.hpp"
#include "../alphabet_container.hpp"
#include "rna4.hpp"

// ------------------------------------------------------------------
// containers
// -----------------------------------------------------------------

namespace seqan3
{

using rna4_vector = std::vector<rna4>;

using rna4_string = std::basic_string<rna4, std::char_traits<rna4>>;

} // namespace seqan3

// ------------------------------------------------------------------
// literals
// -----------------------------------------------------------------

namespace seqan3::literal
{

inline rna4_vector operator "" _rna4(const char * s, std::size_t n)
{
    rna4_vector r;
    r.resize(n);

    std::transform(s, s + n, r.begin(), [] (const char & c)
    {
        return rna4{rna4::char_to_value[c]};
    });

    return r;
}

inline rna4_string operator "" _rna4s(const char * s, std::size_t n)
{
    rna4_string r;
    r.resize(n);

    std::transform(s, s + n, r.begin(), [] (const char & c)
    {
        return rna4{rna4::char_to_value[c]};
    });

    return r;
}

} // namespace seqan3::literal


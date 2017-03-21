#pragma once

#include <iostream>
#include <string>
#include <vector>

#include "../alphabet.hpp"
#include "../alphabet_sequence.hpp"
#include "dna5.hpp"

// ------------------------------------------------------------------
// containers
// -----------------------------------------------------------------

namespace seqan3
{
    
    using dna5_vector = std::vector<dna5>;
    
    using dna5_string = std::basic_string<dna5, std::char_traits<dna5>>;
    
} // namespace seqan3

// ------------------------------------------------------------------
// literals
// -----------------------------------------------------------------

namespace seqan3::literal
{
    
    inline dna5_vector operator "" _dna5(const char * s, std::size_t n)
    {
        dna5_vector r;
        r.resize(n);
        
        std::transform(s, s + n, r.begin(), [] (const char & c)
                       {
                           return dna5{dna5::char_to_value[c]};
                       });
        
        return r;
    }
    
    inline dna5_string operator "" _dna5s(const char * s, std::size_t n)
    {
        dna5_string r;
        r.resize(n);
        
        std::transform(s, s + n, r.begin(), [] (const char & c)
                       {
                           return dna5{dna5::char_to_value[c]};
                       });
        
        return r;
    }
    
} // namespace seqan3::literal


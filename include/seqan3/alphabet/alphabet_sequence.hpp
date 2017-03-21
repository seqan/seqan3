#pragma once

#include <algorithm>
#include <iostream>
#include <string>

#include "alphabet.hpp"

namespace seqan3
{

// ------------------------------------------------------------------
// ostream operator
// ------------------------------------------------------------------

template <typename alphabet_type>
    requires alphabet_concept<alphabet_type>
std::ostream& operator<<(std::ostream & os, std::vector<alphabet_type> const & str)
{
    for (auto c : str)
        os << c;
    return os;
}

template <typename alphabet_type>
    requires alphabet_concept<alphabet_type>
std::ostream& operator<<(std::ostream & os, std::basic_string<alphabet_type> const & str)
{
    for (auto c : str)
        os << c;
    return os;
}


//TODO serialization

}

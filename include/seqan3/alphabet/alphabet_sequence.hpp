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

// ------------------------------------------------------------------
// convert to char container
// ------------------------------------------------------------------

//TODO do we need/want this?
template <typename container_type>
    requires alphabet_concept<typename container_type::value_type>
std::string to_string(container_type const & in)
{
    using alphabet_type = typename container_type::value_type;
    std::string out;
    out.resize(out, in.size());
    std::transform(in.begin(), in.end(), out.begin(), [] (auto c)
    {
        return alphabet_type::value_to_char[static_cast<uint8_t>(c)];
    });
    return out;
}

//TODO serialization

}

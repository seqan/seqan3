#pragma once

#include <iostream>
#include <string>

namespace seqan3
{

// ------------------------------------------------------------------
// concept
// ------------------------------------------------------------------

template <typename t>
concept bool alphabet_concept = requires (t t1, t t2)
{
    // StL concepts
//     requires std::is_pod_v<t> == true;
    requires std::is_swappable_v<t> == true;

    // actual data member (is this required?)
//    t1.value;

    // static data members
    t::value_size;

    // conversion from/to char
    t1.to_char();
    t1.to_integral();

    { t1.from_char(0)     } -> t;
    { t1.from_integral(0) } -> t;


    // required comparison operators
//     { t1 == t2 } -> bool;
//     { t1 != t2 } -> bool;
//     { t1 <  t2 } -> bool;
//     { t1 >  t2 } -> bool;
//     { t1 <= t2 } -> bool;
//     { t1 >= t2 } -> bool;
};

// ------------------------------------------------------------------
// ostream operator
// ------------------------------------------------------------------

template <typename alphabet_type>
    requires alphabet_concept<alphabet_type>
std::ostream& operator<<(std::ostream & os, alphabet_type const & c)
{
    os << c.to_char();
    return os;
}

//TODO serialization

} // namespace seqan3


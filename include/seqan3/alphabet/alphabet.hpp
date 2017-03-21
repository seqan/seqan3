#pragma once

#include <iostream>
#include <string>

namespace seqan3
{

// ------------------------------------------------------------------
// free functions to member function wrapper
// ------------------------------------------------------------------
// move to alphabet_detail.hpp?


/* The public alphabet concept requires only free function access
 * for type that have member functions we create a wrapper here
 * so you don't have to repeat it.
 */

namespace detail
{

template <typename t>
concept bool internal_alphabet_concept = requires (t t1)
{
    typename t::char_type;
    typename t::integral_type;
    t::value_size;

    { t1.to_char() } -> typename t::char_type;
    { t1.to_integral() } -> typename t::integral_type;

    { t1.from_char('a')   } -> t;
    { t1.from_integral(0) } -> t;
};
}

template <typename alphabet_type>
     requires detail::internal_alphabet_concept<alphabet_type>
constexpr auto value_size(alphabet_type const &)
{
    return alphabet_type::value_size;
}

template <typename alphabet_type>
    requires detail::internal_alphabet_concept<alphabet_type>
constexpr char to_char(alphabet_type const & c)
{
    return c.to_char();
}

template <typename alphabet_type>
    requires detail::internal_alphabet_concept<alphabet_type>
constexpr auto to_integral(alphabet_type const & c)
{
    return c.to_integral();
}

template <typename alphabet_type>
    requires detail::internal_alphabet_concept<alphabet_type>
constexpr alphabet_type from_char(alphabet_type & c, char const in)
{
    return c.from_char(in);
}

template <typename alphabet_type, typename input_type>
    requires detail::internal_alphabet_concept<alphabet_type> &&
             std::is_integral<input_type>::value
constexpr alphabet_type from_integral(alphabet_type & c, input_type const in)
{
    return c.from_integral(in);
}

template <typename alphabet_type>
    requires detail::internal_alphabet_concept<alphabet_type>
struct alphabet_char
{
    using type = typename alphabet_type::char_type;
};

template <typename alphabet_type>
//     requires detail::internal_alphabet_concept<alphabet_type>
using alphabet_char_t = typename alphabet_char<alphabet_type>::type;

template <typename alphabet_type>
    requires detail::internal_alphabet_concept<alphabet_type>
struct alphabet_integral
{
    using type = typename alphabet_type::integral_type;
};

template <typename alphabet_type>
//     requires detail::internal_alphabet_concept<alphabet_type>
using alphabet_integral_t = typename alphabet_integral<alphabet_type>::type;

// ------------------------------------------------------------------
// alphabet concept
// ------------------------------------------------------------------

template <typename t>
concept bool alphabet_concept = requires (t t1, t t2)
{
    // StL concepts
    requires std::is_pod_v<t> == true;
    requires std::is_swappable_v<t> == true;

    // static data members
    { value_size(t1) };

    // conversion from/to char
    { to_char(t1)     } -> alphabet_char_t<t>;
    { to_integral(t1) } -> alphabet_integral_t<t>;

    { from_char(t1, 0)     } -> t;
    { from_integral(t1, 0) } -> t;

    // required comparison operators
    { t1 == t2 } -> bool;
    { t1 != t2 } -> bool;
    { t1 <  t2 } -> bool;
    { t1 >  t2 } -> bool;
    { t1 <= t2 } -> bool;
    { t1 >= t2 } -> bool;
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

#include <iostream>

//! [usage]
#include <seqan3/core/detail/reflection.hpp>
#include <seqan3/core/debug_stream.hpp>

namespace foo
{

template <typename ...type>
struct bar
{};

} // namespace foo

int main()
{
    seqan3::debug_stream << seqan3::detail::get_display_name_v<foo::bar<char, double>> << std::endl; // prints: foo::bar<char, double> >
    return 0;
}
//! [usage]

#include <iostream>

//! [usage]
#include <seqan3/core/detail/reflection.hpp>

namespace foo
{

template <typename ...type>
struct bar
{};

} // namespace foo

int main()
{
    std::cout << seqan3::detail::get_display_name_v<foo::bar<char, double>>.string() << std::endl; // prints: foo::bar<char, double> >
    return 0;
}
//! [usage]

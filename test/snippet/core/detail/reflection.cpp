#include <iostream>

#include <seqan3/core/detail/reflection.hpp>

using namespace seqan3;

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
    std::cout << detail::get_display_name_v<foo::bar<char, double>>.string() << std::endl; // prints: foo::bar<char, double> >
}
//! [usage]

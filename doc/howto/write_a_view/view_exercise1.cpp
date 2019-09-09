//![start]
#include <iostream>
#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/std/ranges>

using namespace seqan3;

//![start]
auto const my_convert_to_char_view = std::views::transform([] (auto const alph)
{
    return to_char(alph);
});

//![end]
int main()
{
    std::vector<dna5> vec{"ATTAGATTA"_dna5};
    // std::cout << vec[0] << '\n';                 // won't work

    auto v = vec | my_convert_to_char_view;

    std::cout << v[0] << '\n';                      // prints "A"
}
//![end]

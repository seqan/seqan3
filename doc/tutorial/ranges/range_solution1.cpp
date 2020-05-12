#include <iostream>
#include <vector>
#include <seqan3/std/ranges>                    // include all of the standard library's views

int main()
{
    std::vector vec{1, 2, 3, 4, 5, 6};
    auto v = vec | std::views::filter([] (auto const i) { return i % 2 == 0; })
                 | std::views::transform([] (auto const i) { return i*i; });

    std::cout << *v.begin() << '\n'; // prints 4
}

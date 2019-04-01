#include <iostream>
#include <vector>
#include <seqan3/std/ranges> // include all of the standard library's views

int main()
{
{
//![all]
//![def]
std::vector vec{1, 2, 3, 4, 5, 6};
//![rev_def]
auto v = std::view::reverse(vec);
//![rev_def]
//![def]

std::cout << *v.begin() << '\n';
//![all]
}

{
//![assign_through]
//![piped]
std::vector vec{1, 2, 3, 4, 5, 6};
auto v = vec | std::view::reverse | std::view::drop(2);

std::cout << *v.begin() << '\n';
//![piped]
*v.begin() = 42;                  // now vec == {1, 2, 3, 42, 5, 6 } !!
//![assign_through]
}

{
//![solution1]
std::vector vec{1, 2, 3, 4, 5, 6};
auto v = vec | std::view::filter(   [] (auto const i) { return i % 2 == 0; })
             | std::view::transform([] (auto const i) { return i*i; });

std::cout << *v.begin() << '\n'; // prints 4
//![solution1]
}
}

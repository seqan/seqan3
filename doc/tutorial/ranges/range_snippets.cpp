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
auto v = std::views::reverse(vec);
//![rev_def]
//![def]

std::cout << *v.begin() << '\n';
//![all]
}

{
//![assign_through]
//![piped]
std::vector vec{1, 2, 3, 4, 5, 6};
auto v = vec | std::views::reverse | std::views::drop(2);

std::cout << *v.begin() << '\n';
//![piped]
*v.begin() = 42;                  // now vec == {1, 2, 3, 42, 5, 6 } !!
//![assign_through]
}
}

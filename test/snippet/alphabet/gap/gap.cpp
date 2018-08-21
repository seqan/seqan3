#include <seqan3/alphabet/gap/gap.hpp>

using namespace seqan3;

int main()
{
//! [general]
gap my_gap = gap::GAP;
gap another_gap{};
another_gap.assign_char('A'); // setting this does not change anything

std::cout << my_gap.to_char(); // outputs '-'
if (my_gap.to_char() == another_gap.to_char())
    std::cout << "Both gaps are the same!";
//! [general]
}

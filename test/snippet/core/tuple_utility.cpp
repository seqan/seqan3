#include <seqan3/core/tuple_utility.hpp>

using namespace seqan3;

int main()
{
//! [usage]
// Split at position 2.
std::tuple<int, char, float, std::string> t{1, 'c', 0.3, "hello"};
auto [left, right] = tuple_split<2>(t);  // decltype(left) -> std::tuple<int, char>; decltype(right) -> std::tuple<float, std::string>;

// Split at position 0.
auto [left1, right1] = tuple_split<0>(t);  // decltype(left1) -> std::tuple<>; decltype(right1) -> std::tuple<int, char, float, std::string>;

// Split at position 4.
auto [left2, right2] = tuple_split<4>(t);  // decltype(left1) -> std::tuple<int, char, float, std::string>; decltype(right1) -> std::tuple<>;
//! [usage]
(void) left;
(void) left1;
(void) left2;
(void) right;
(void) right1;
(void) right2;
}

#include <seqan3/utility/tuple/split.hpp>

int main()
{
    // Split at position 2.
    std::tuple<int, char, float, std::string> t{1, 'c', 0.3, "hello"};
    auto [left, right] = seqan3::tuple_split<2>(t);
    // decltype(left) -> std::tuple<int, char>; decltype(right) -> std::tuple<float, std::string>;

    // Split at position 0.
    auto [left1, right1] = seqan3::tuple_split<0>(t);
    // decltype(left1) -> std::tuple<>; decltype(right1) -> std::tuple<int, char, float, std::string>;

    // Split at position 4.
    auto [left2, right2] = seqan3::tuple_split<4>(t);
    // decltype(left2) -> std::tuple<int, char, float, std::string>; decltype(right2) -> std::tuple<>;
}

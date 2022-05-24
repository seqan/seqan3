#include <iostream>

#include <seqan3/utility/char_operations/predicate.hpp>

int main()
{
    char chr{'1'};
    constexpr auto my_cond = seqan3::is_char<'%'> || seqan3::is_digit;
    bool is_percent = my_cond(chr);
    std::cout << std::boolalpha << is_percent << '\n'; // true
}

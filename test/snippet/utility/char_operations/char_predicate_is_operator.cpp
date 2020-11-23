#include <seqan3/utility/char_operations/predicate.hpp>

int main()
{
    char chr{'1'};
    auto constexpr my_cond = seqan3::is_char<'%'> || seqan3::is_digit;
    bool is_percent = my_cond(chr); // is_percent == true
}

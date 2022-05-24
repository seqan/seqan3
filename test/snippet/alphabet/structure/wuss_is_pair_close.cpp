#include <iostream>

#include <seqan3/alphabet/structure/wuss.hpp>

int main()
{
    using namespace seqan3::literals;

    bool is_closing_char_member = '}'_wuss51.is_pair_close();
    bool is_closing_char_free = seqan3::is_pair_close('.'_wuss51);

    std::cout << std::boolalpha << is_closing_char_member << '\n'; // true
    std::cout << std::boolalpha << is_closing_char_free << '\n';   // false
}

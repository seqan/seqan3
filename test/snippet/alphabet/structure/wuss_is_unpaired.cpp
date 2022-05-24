#include <iostream>

#include <seqan3/alphabet/structure/wuss.hpp>

int main()
{
    using namespace seqan3::literals;

    bool is_unpaired_char_member = '.'_wuss51.is_unpaired();
    bool is_unpaired_char_free = seqan3::is_unpaired('{'_wuss51);

    std::cout << std::boolalpha << is_unpaired_char_member << '\n'; // true
    std::cout << std::boolalpha << is_unpaired_char_free << '\n';   // false
}

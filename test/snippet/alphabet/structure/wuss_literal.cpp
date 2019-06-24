#include <seqan3/alphabet/structure/wuss.hpp>

int main()
{
    using seqan3::operator""_wuss51;

    // Using the string literal to assign a vector of WUSS characters:
    std::vector<seqan3::wuss<>> foo{".<..>."_wuss51};
    std::vector<seqan3::wuss<>> bar = ".<..>."_wuss51;
    auto bax = ".<..>."_wuss51;

    // Using the char literal to assign a single WUSS character:
    seqan3::wuss51 my_letter{'~'_wuss51};
    // does not work:
    // wuss51 my_letter{'~'}; // <- char not implicitly convertible
    // works for each wuss alphabet size:
    my_letter.assign_char('<'); // <- assigns the char explicitly
}

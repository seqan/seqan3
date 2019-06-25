#include <seqan3/alphabet/structure/wuss.hpp>

int main()
{
    using seqan3::operator""_wuss51;

    // Using the char literal to assign a single WUSS character:
    seqan3::wuss51 my_letter{'~'_wuss51};
    
    my_letter.assign_char('<'); // <- assigns the char explicitly
}

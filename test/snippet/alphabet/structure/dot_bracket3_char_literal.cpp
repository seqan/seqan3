#include <seqan3/alphabet/structure/dot_bracket3.hpp>

int main()
{
    using seqan3::operator""_db3;

    // Using the char literal to assign a single dot bracket:
    seqan3::dot_bracket3 my_letter{'('_db3};
    
    my_letter.assign_char(')');             // <- assigns the char explicitly
}

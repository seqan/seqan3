#include <seqan3/alphabet/structure/dot_bracket3.hpp>

int main()
{
    using seqan3::operator""_db3;
    
    // Using the string literal to assign a vector of dot brackets:
    std::vector<seqan3::dot_bracket3> foo{".(..)."_db3};
    std::vector<seqan3::dot_bracket3> bar = ".(..)."_db3;
    auto bax = ".(..)."_db3;

    // Using the char literal to assign a single dot bracket:
    seqan3::dot_bracket3 my_letter{'('_db3};
    // does not work:
    // seqan3::dot_bracket3 my_letter{'('}; // <- char not implicitly convertible
    my_letter.assign_char(')');     // <- assigns the char explicitly
}

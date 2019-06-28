#include <vector>

#include <seqan3/alphabet/structure/dot_bracket3.hpp>

int main()
{
    using seqan3::operator""_db3;

    // Using the string literal to assign a vector of dot brackets:
    std::vector<seqan3::dot_bracket3> foo{".(..)."_db3};
    std::vector<seqan3::dot_bracket3> bar = ".(..)."_db3;
    auto bax = ".(..)."_db3;
}

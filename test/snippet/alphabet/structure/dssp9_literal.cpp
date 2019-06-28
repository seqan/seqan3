#include <vector>

#include <seqan3/alphabet/structure/dssp9.hpp>

int main()
{
    using seqan3::operator""_dssp9;

    // Using the string literal to assign a vector of DSSP annotations:
    std::vector<seqan3::dssp9> foo{"EHHHHT"_dssp9};
    std::vector<seqan3::dssp9> bar = "EHHHHT"_dssp9;
    auto bax = "EHHHHT"_dssp9;
}

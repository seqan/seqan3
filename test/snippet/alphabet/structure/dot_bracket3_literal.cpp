#include <seqan3/alphabet/structure/dot_bracket3.hpp>

int main()
{
    using namespace seqan3::literals;

    std::vector<seqan3::dot_bracket3> sequence1{".(..)."_db3};
    std::vector<seqan3::dot_bracket3> sequence2 = ".(..)."_db3;
    auto sequence3 = ".(..)."_db3;
}

#include <seqan3/alphabet/structure/dot_bracket3.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::dot_bracket3 letter1{'('_db3};
    auto letter2 = '('_db3;
}

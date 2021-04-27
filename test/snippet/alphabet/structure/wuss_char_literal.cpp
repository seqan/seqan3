#include <seqan3/alphabet/structure/wuss.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::wuss51 letter1{'('_wuss51};
    auto letter2 = '('_wuss51;
}

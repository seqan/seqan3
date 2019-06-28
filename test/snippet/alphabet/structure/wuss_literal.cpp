#include <seqan3/alphabet/structure/wuss.hpp>

int main()
{
    using seqan3::operator""_wuss51;

    // Using the string literal to assign a vector of WUSS characters:
    std::vector<seqan3::wuss<>> foo{".<..>."_wuss51};
    std::vector<seqan3::wuss<>> bar = ".<..>."_wuss51;
    auto bax = ".<..>."_wuss51;
}

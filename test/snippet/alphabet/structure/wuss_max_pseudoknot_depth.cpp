#include <iostream>

#include <seqan3/alphabet/structure/wuss.hpp>

int main()
{
    uint8_t max_depth_member = seqan3::wuss51::max_pseudoknot_depth;
    uint8_t max_depth_meta = seqan3::max_pseudoknot_depth<seqan3::wuss51>;
    std::cout << static_cast<uint16_t>(max_depth_member) << '\n'; // 22
    std::cout << static_cast<uint16_t>(max_depth_meta) << '\n';   // 22
}

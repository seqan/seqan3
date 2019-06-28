#include <seqan3/alphabet/structure/wuss.hpp>

int main()
{
    uint8_t max_depth = seqan3::wuss51::max_pseudoknot_depth;             // 22
    uint8_t max_depth_alt = seqan3::max_pseudoknot_depth<seqan3::wuss51>; // 22
}

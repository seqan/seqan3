#include <seqan3/alphabet/structure/wuss.hpp>

int main()
{
    using namespace seqan3::literals;

    bool is_unpaired_char = '.'_wuss51.is_unpaired();                // true
    bool is_unpaired_char_alt = seqan3::is_unpaired('{'_wuss51);     // false
}

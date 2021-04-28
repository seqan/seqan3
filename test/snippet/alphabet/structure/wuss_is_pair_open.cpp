#include <seqan3/alphabet/structure/wuss.hpp>

int main()
{
    using namespace seqan3::literals;

    bool is_opening_char = '{'_wuss51.is_pair_open();                // true
    bool is_opening_char_alt = seqan3::is_pair_open('.'_wuss51);     // false
}

#include <seqan3/alphabet/structure/wuss.hpp>

int main()
{
    using namespace seqan3::literals;

    bool is_closing_char = '}'_wuss51.is_pair_close();               // true
    bool is_closing_char_alt = seqan3::is_pair_close('.'_wuss51);    // false
}

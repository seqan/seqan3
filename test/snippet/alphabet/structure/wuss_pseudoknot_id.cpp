#include <seqan3/alphabet/structure/wuss.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using namespace seqan3::literals;

    auto pk_opt = '.'_wuss51.pseudoknot_id();                       // std::optional -> false
    pk_opt = seqan3::pseudoknot_id('{'_wuss51);                     // std::optional -> true: 3

    if (pk_opt)
        seqan3::debug_stream << *pk_opt << '\n';                    // 3
}

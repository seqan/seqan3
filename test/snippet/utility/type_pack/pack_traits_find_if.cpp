#include <seqan3/utility/type_pack/traits.hpp>

int main()
{
    // None of the types in t is a pointer so find_if returns -1. However, int and bool are both integral,
    // so find_if returns 0 for the first occurrence.
    static_assert(seqan3::pack_traits::find_if<std::is_pointer, int, float, double> == -1);
    static_assert(seqan3::pack_traits::find_if<std::is_integral, int, float, double> == 0);
}

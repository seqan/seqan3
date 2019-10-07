#include <seqan3/core/type_list/traits.hpp>

int main()
{
    // Get the size of the pack.
    static_assert(seqan3::pack_traits::size<int, float, bool, int> == 4);
}

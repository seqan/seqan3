#include <seqan3/utility/type_pack/traits.hpp>

int main()
{
    // Get the size of the pack.
    static_assert(seqan3::pack_traits::size<int, float, bool, int> == 4);
}

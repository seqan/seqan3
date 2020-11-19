#include <seqan3/utility/type_pack/traits.hpp>

int main()
{
    // Count the number of type int in the pack.
    static_assert(seqan3::pack_traits::count<int, int, float, bool, int> == 2);
}

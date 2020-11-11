#include <seqan3/utility/type_pack/traits.hpp>

int main()
{
    // Take the last two types in the pack.
    static_assert(std::same_as<seqan3::type_list<bool, int>, seqan3::pack_traits::take_last<2, int, float, bool, int>>);
}

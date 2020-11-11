#include <seqan3/utility/type_pack/traits.hpp>

int main()
{
    // Drop the last two types in the pack.
    static_assert(std::same_as<seqan3::type_list<int, float>, seqan3::pack_traits::drop_last<2, int, float, bool, int>>);
}

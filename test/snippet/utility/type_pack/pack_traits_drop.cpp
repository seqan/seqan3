#include <seqan3/utility/type_pack/traits.hpp>

int main()
{
    // Drop the first two types in the pack.
    static_assert(std::same_as<seqan3::type_list<bool, int>, seqan3::pack_traits::drop<2, float, double, bool, int>>);
}

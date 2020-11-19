#include <seqan3/utility/type_pack/traits.hpp>

int main()
{
    // Check if the first value is int.
    static_assert(std::same_as<int, seqan3::pack_traits::front<int, float, bool, int, float>>);
}

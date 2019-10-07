#include <seqan3/core/type_list/traits.hpp>

int main()
{
    // Check if the first value is int.
    static_assert(std::same_as<int, seqan3::pack_traits::front<int, float, bool, int, float>>);
}

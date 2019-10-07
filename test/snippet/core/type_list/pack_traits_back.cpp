#include <seqan3/core/type_list/traits.hpp>

int main()
{
    // Check if the last value is float.
    static_assert(std::same_as<float, seqan3::pack_traits::back<int, float, bool, int, float>>);
}

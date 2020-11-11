#include <seqan3/utility/type_pack/traits.hpp>

int main()
{
    // Check if the last value is float.
    static_assert(std::same_as<float, seqan3::pack_traits::back<int, float, bool, int, float>>);
}

#include <seqan3/utility/type_list/traits.hpp>

int main()
{
    using list_t = seqan3::type_list<int, float, bool, int, float>;

    // Access the last value (float) with seqan3::list_traits::back
    static_assert(std::same_as<float, seqan3::list_traits::back<list_t>>);
}

#include <seqan3/utility/type_list/traits.hpp>

int main()
{
    using list_t = seqan3::type_list<int, float, bool, int, float>;

    static_assert(std::same_as<int, seqan3::list_traits::front<list_t>>); // Check if the first value is int.
}

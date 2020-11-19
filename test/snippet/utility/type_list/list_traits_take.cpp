#include <seqan3/utility/type_list/traits.hpp>

int main()
{
    using list_t = seqan3::type_list<int, float, bool, int>;

    // Take the first two types in list_t.
    static_assert(std::same_as<seqan3::type_list<int, float>, seqan3::list_traits::take<2, list_t>>);
}

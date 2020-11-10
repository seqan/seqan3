#include <seqan3/utility/type_list/traits.hpp>

int main()
{
    using list_t = seqan3::type_list<int, float, bool, int>;

    // Drop the first two types in list_t.
    static_assert(std::same_as<seqan3::type_list<bool, int>, seqan3::list_traits::drop<2, list_t>>);
}

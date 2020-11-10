#include <seqan3/utility/type_list/traits.hpp>

int main()
{
    using list_t = seqan3::type_list<int, float, bool, int>;

    // Drop the last two types in list_t.
    static_assert(std::same_as<seqan3::type_list<int, float>, seqan3::list_traits::drop_last<2, list_t>>);
}

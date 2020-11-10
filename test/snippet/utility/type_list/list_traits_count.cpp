#include <seqan3/utility/type_list/traits.hpp>

int main()
{
    using list_t = seqan3::type_list<int, float, bool, int>;

    // Count the number of type int in list_t.
    static_assert(seqan3::list_traits::count<int, list_t> == 2);
}

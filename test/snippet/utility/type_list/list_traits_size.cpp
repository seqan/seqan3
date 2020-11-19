#include <seqan3/utility/type_list/traits.hpp>

int main()
{
    using list_t = seqan3::type_list<int, float, bool, int>;

    static_assert(seqan3::list_traits::size<list_t> == 4);
}

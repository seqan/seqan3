#include <seqan3/utility/type_list/traits.hpp>

int main()
{
    using list_t = seqan3::type_list<int, float, bool, int>;

    static_assert(std::same_as<seqan3::type_list<float, bool, int>, seqan3::list_traits::drop_front<list_t>>);
}

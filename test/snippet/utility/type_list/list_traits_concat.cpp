#include <seqan3/utility/type_list/traits.hpp>

int main()
{
    using list_t = seqan3::type_list<int, float, bool>;
    using list_t2 = seqan3::type_list<double, char, int>;
    using list_t3 = seqan3::type_list<int, int>;

    static_assert(std::same_as<seqan3::list_traits::concat<list_t, list_t2, list_t3>,
                               seqan3::type_list<int, float, bool, double, char, int, int, int>>);
}

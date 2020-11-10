#include <seqan3/utility/type_list/traits.hpp>

int main()
{
    using list_t = seqan3::type_list<int, float, bool, double, char, int>;
    using split_t = seqan3::list_traits::split_after<3, list_t>;

    // Use ::first_type and ::second_type to access the type lists after being split.
    static_assert(std::same_as<seqan3::type_list<int, float, bool>,
                               split_t::first_type>);
    static_assert(std::same_as<seqan3::type_list<double, char, int>,
                               split_t::second_type>);
}

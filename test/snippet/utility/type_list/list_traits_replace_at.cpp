#include <seqan3/utility/type_list/traits.hpp>

int main()
{
    using list_t = seqan3::type_list<int, float, bool, double>;

    // Replace the second element with int.
    static_assert(std::same_as<seqan3::type_list<int, int, bool, double>, seqan3::list_traits::replace_at<int, 1, list_t>>);
}

#include <seqan3/utility/type_list/traits.hpp>

int main()
{
    using list_t = seqan3::type_list<int, float, bool, double>;

    // Look at the 2nd element.
    static_assert(std::same_as<float, seqan3::list_traits::at<1, list_t>>);
    // Look at the last element.
    static_assert(std::same_as<double, seqan3::list_traits::at<-1, list_t>>);
}

#include <seqan3/utility/type_pack/traits.hpp>

int main()
{
    // Look at the 2nd element.
    static_assert(std::same_as<float, seqan3::pack_traits::at<1, int, float, bool, double>>);
    // Look at the last element.
    static_assert(std::same_as<double, seqan3::pack_traits::at<-1, int, float, bool, double>>);
}

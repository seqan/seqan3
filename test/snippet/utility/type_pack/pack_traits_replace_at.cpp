#include <seqan3/utility/type_pack/traits.hpp>

int main()
{
    // Replace the second element of the type pack with int.
    static_assert(std::same_as<seqan3::type_list<int, int, bool, double>,
                               seqan3::pack_traits::replace_at<int, 1, int, float, bool, double>>);
}

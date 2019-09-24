#include <list>
#include <vector>

#include <seqan3/core/type_list/traits.hpp>
#include <seqan3/core/type_traits/range.hpp>

int main()
{
    // Transform the types in the pack into reference types.
    static_assert(std::same_as<seqan3::pack_traits::transform<seqan3::reference_t, std::vector<int>,
                                                                                   std::vector<float>,
                                                                                   std::list<bool>>,
                               seqan3::type_list<int &, float &, bool &>>);
}

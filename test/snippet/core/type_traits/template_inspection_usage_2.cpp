#include <vector>

#include <seqan3/core/type_traits/template_inspection.hpp>
#include <seqan3/core/type_list/type_list.hpp>

int main()
{
    using my_type = std::vector<int>;

    if constexpr (seqan3::detail::is_type_specialisation_of_v<my_type, std::vector>) // std::vector has no <> !
    {
        // ...
    }
}

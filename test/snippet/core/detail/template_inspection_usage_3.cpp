#include <vector>

#include <seqan3/core/detail/template_inspection.hpp>

int main()
{
    using my_type = std::vector<int>;

    if constexpr (seqan3::detail::template_specialisation_of<my_type, std::vector>) // Note: std::vector has no <> !
    {
        // ...
    }
}

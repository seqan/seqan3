#include <seqan3/core/type_traits/template_inspection.hpp>
#include <seqan3/core/algorithm/configuration.hpp>

using namespace seqan3;

int main()
{
//! [usage]
using tl = type_list<int, char, double>;
using t = detail::transfer_template_args_onto_t<tl, std::tuple>;
// t is std::tuple<int, char, double>
//! [usage]

// suppresses the warning that the type t is unused
[[maybe_unused]] t test;

//! [usage_2]
using my_type = std::vector<int>;

if constexpr (detail::is_type_specialisation_of_v<my_type, std::vector>) // std::vector has no <> !
{
    // ...
}
//! [usage_2]
}

#include <seqan3/core/type_traits/template_inspection.hpp>
#include <seqan3/core/type_list/type_list.hpp>

int main()
{
    using tl = seqan3::type_list<int, char, double>;
    using t = seqan3::detail::transfer_template_args_onto_t<tl, std::tuple>;
    // t is std::tuple<int, char, double>
}

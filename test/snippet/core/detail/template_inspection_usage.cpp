#include <seqan3/core/detail/template_inspection.hpp>
#include <seqan3/utility/type_list/type_list.hpp>

int main()
{
    using list_to_transfer = seqan3::type_list<int, char, double>;
    using resulting_t = seqan3::detail::transfer_template_args_onto_t<list_to_transfer, std::tuple>;

    static_assert(std::same_as<resulting_t, std::tuple<int, char, double>>);
}

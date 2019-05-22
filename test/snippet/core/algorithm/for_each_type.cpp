#include <iostream>
#include <string>

#include <seqan3/core/algorithm/parameter_pack.hpp>
#include <seqan3/core/debug_stream.hpp>

using namespace seqan3;

namespace incomplete
{

struct type;

} // namespace incomplete

int main()
{
    // With c++20 you could also write it like this
    // auto fn = []<typename type>(std::type_identity<type>)
    // {
    // ...
    // };
    auto fn = [](auto id)
    {
        // id is of type std::type_identity<type>
        using id_t = decltype(id);
        using type = typename id_t::type;

        static_assert(std::is_same_v<id_t, std::type_identity<type>>, "id is of type std::type_identity<type>");

        if constexpr(std::is_same_v<type, bool>)
            debug_stream << "bool";
        else if constexpr(std::is_same_v<type, int>)
            debug_stream << "int";
        else if constexpr(std::is_same_v<type, float>)
            debug_stream << "float";
        else if constexpr(std::is_same_v<type, incomplete::type>)
            debug_stream << "incomplete::type";

        debug_stream << ", ";
    };

    // prints each type name, i.e. "int, float, bool, incomplete::type, \n"
    detail::for_each_type<int, float, bool, incomplete::type>(fn);
    debug_stream << "\n";

    // is the same as explicitly writing
    fn(std::type_identity<int>{});
    fn(std::type_identity<float>{});
    fn(std::type_identity<bool>{});
    fn(std::type_identity<incomplete::type>{});
    debug_stream << "\n";
    return 0;
}

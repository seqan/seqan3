//! [all]
#include <variant> // for std::visit

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/alignment_file/sam_tag_dictionary.hpp>
#include <seqan3/range/container/concept.hpp> // for the seqan3::Container

using namespace seqan3;

// a lambda helper function that prints every type in the std::variant<...allowed SAM tag types...>
auto print_fn = [] (auto && arg)
{
    using T = remove_cvref_t<decltype(arg)>; // the type T of arg.

    if constexpr (!Container<T>)     // If T is not a container,
    {
        debug_stream << arg << std::endl;       // just print arg directly.
    }
    else                                     // If T is a container,
    {
        for (auto const & arg_v : arg)       // print every value in arg.
            debug_stream << arg_v << ",";
        debug_stream << std::endl;
    }
};

int main()
{
    sam_tag_dictionary dict{};          // initialise empty dictionary

    // ! there is no get function for unknown tags !
    // dict.get<"XZ"_tag>() = 3;
    // but you can use the operator[]
    dict["XZ"_tag] = 3; // set unknown SAM tag 'XZ' to 3 (type int32_t)

    // ! there is no get function for unknown tags !
    // auto nm = dict.get<"XZ"_tag>();
    // but you can use the operator[] again
    auto xz = dict["XZ"_tag]; // get SAM tag 'XZ' (type std::variant<...allowed SAM tag types...>)

    // ! you cannot print a std::variant directly !
    // debug_stream << nm << std::endl;
    // but you can use visit:
    std::visit(print_fn, xz); // prints 3
}
//! [all]

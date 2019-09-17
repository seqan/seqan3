#include <type_traits>

#include <seqan3/core/type_traits/lazy.hpp>
#include <seqan3/core/debug_stream.hpp>

template <typename t>
    requires std::is_integral_v<t>
struct foo
{
    t value;
};

template <typename t>
auto bar(t const & v)
{
    using cond_t = seqan3::detail::lazy_conditional_t<seqan3::detail::is_instantiable_with_v<foo, t>,
                                                      seqan3::detail::lazy<foo, t>,
                                                      t>;
    return cond_t{v};
}

int main()
{
    foo<int> a = bar(10);
    seqan3::debug_stream << "a: " << a.value << '\n'; // prints 10
    float b = bar(0.4f);
    seqan3::debug_stream << "b: " << b << '\n'; // prints 0.4
}

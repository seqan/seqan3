#include <seqan3/range/view/detail.hpp>

using namespace seqan3;
using namespace seqan3::detail;

template <typename urng_t>
requires std::ranges::InputRange<urng_t> // or more/stricter requirements
struct view_foo
{
    urng_t const & urange;
    int param;

    view_foo(urng_t const & _urange) :
        urange{_urange}
    {}

    view_foo(urng_t const & _urange, int const & _param) :
        urange{_urange}, param {_param}
    {}

    view_foo(urng_t const & _urange, int && _param) :
        urange{_urange}, param{std::move(_param)}
    {}
};

using foo_fn = generic_pipable_view_adaptor<view_foo>;

namespace view
{
    inline constexpr foo_fn foo;
}
#if 0 // Copied code for documentation.
//! [usage]
template <typename urng_t>
    requires std::ranges::InputRange<urng_t> // or more/stricter requirements
struct view_foo
{
    // your implementation
};

using foo_fn = generic_pipable_view_adaptor<view_foo>;

namespace view
{
    inline constexpr foo_fn foo;
}
//! [usage]
#endif
int main()
{
{
//! [function_call]
std::vector<int> container{1, 2, 3};

auto w = view::foo(container);    // if the view takes no constructor args beyond urange
auto v = view::foo(container, 7); // if the view takes e.g. an extra int argument
// in both cases v is now of type view_foo<std::vector<int>>
//! [function_call]
(void) v;
(void) w;
}

{
//! [function_call_2]
std::vector<int> container{1, 2, 3};

auto v = view::foo(7); // v is NOT OF TYPE view_foo<std::vector<int>>

// it's usually not used like above, instead use it inside a pipe:
auto w = container | view::foo(7);
//! [function_call_2]
(void) v;
(void) w;
}

{
//! [pipe_syntax]
std::vector<int> container{1, 2, 3};

auto v = container | view::foo(7);
//                 ^           ^
//     this operator           the intermediate operator() that returns a bound functor
//! [pipe_syntax]
(void) v;
}

{
//! [pipe_syntax_2]
std::vector<int> container{1, 2, 3};

auto v = container | view::foo; // v is now of type view_foo<std::vector<int>
//! [pipe_syntax_2]
(void) v;
}

#if 0 // Copied code for documentation.
//! [adaptor_template]
// 1.
template <typename urng_t>
    requires std::ranges::InputRange<urng_t> // or more/stricter requirements
struct view_foo
{
// your implementation
};

// 2. this is sufficient to declare the adaptor type:
using foo_fn = detail::generic_pipable_view_adaptor<view_foo>;

namespace view
{
// 3. the adaptor instance
inline constexpr foo_fn foo;
}
//! [adaptor_template]
#endif
}

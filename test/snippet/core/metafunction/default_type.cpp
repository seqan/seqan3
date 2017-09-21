#include <seqan3/core/metafunction/default_type.hpp>

using namespace seqan3::detail;

template <typename T>
struct A;

template <>
struct A<int>
{
    using type = int;
};

// A<unsigned>::type is not defined, thus falling back to `void`
static_assert(std::is_same_v<void, default_type_t<A<unsigned>, void>>);

// A<int>::type is defined, use A<int>::type
static_assert(std::is_same_v<int, default_type_t<A<int>, void>>);

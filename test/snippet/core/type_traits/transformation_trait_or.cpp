#include <seqan3/core/type_traits/transformation_trait_or.hpp>

using namespace seqan3::detail;

template <typename T>
struct A;

template <>
struct A<int>
{
    using type = int;
};

// A<unsigned>::type is not defined, thus falling back to `void`
static_assert(std::is_same_v<void, transformation_trait_or_t<A<unsigned>, void>>);

// A<int>::type is defined, use A<int>::type
static_assert(std::is_same_v<int, transformation_trait_or_t<A<int>, void>>);

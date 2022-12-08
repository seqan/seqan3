#include <cmath>

#include <seqan3/core/debug_stream.hpp>

//![concept]
// helper concept has_foo:
template <typename T>
concept has_foo = requires (T val) {
                      typename T::FOO; // requirement 1
                      val.foo;         // requirement 2
                  };

// concept fooger:
template <typename T>
concept fooger = has_foo<T> && std::same_as<typename T::FOO, int>;
//![concept]

struct my_type
{
    using FOO = int;
    char foo{}; // foo can be of any type, here it is of type `char`
};

//![main]
int main()
{
    seqan3::debug_stream << fooger<my_type> << std::endl; // should print 1
}
//![main]

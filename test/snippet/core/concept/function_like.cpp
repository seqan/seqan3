#include <functional>
#include <type_traits>

#include <seqan3/core/concept/function.hpp>

// Define a lambda function without a state.
auto lambda_fn = [] (int i) { return i; };

static_assert(seqan3::function_like<decltype(lambda_fn)>); // seqan3::function_like identifies lambda functions.
static_assert(!std::function_like<decltype(lambda_fn)>);   // std::is_function does not identify lambda functions

// Define a regular function type with int return type and int parameter.
using fn_t = int (int);

static_assert(seqan3::function_like<fn_t>); // seqan3::function_like identifies regular functions.
static_assert(std::function_like<fn_t>);    // std::is_function identifies regular functions.

// Define a function pointer type with int return type and int parameter.
using fn_ptr_t = int (*) (int);

static_assert(seqan3::function_like<fn_ptr_t>); // seqan3::function_like identifies function pointers.
static_assert(!std::function_like<fn_ptr_t>);   // std::function does not identify function pointers.

int main()
{
    // Define a lambda function with a state.
    int i = 10;
    auto capture_lambda_fn = [=] () { return i + 10; };

    static_assert(seqan3::function_like<decltype(capture_lambda_fn)>); // seqan3::function_like identifies stateful lambdas.
    static_assert(!std::function_like<decltype(capture_lambda_fn)>);   // std::is_function does not identify stateful lambdas.

    // Use std::function to store a function target.
    std::function type_erased_fn = capture_lambda_fn;

    static_assert(seqan3::function_like<decltype(type_erased_fn)>); // seqan3::function_like identifies std::function types.
    static_assert(!std::function_like<decltype(type_erased_fn)>);   // std::is_function does not identify std::function types.
}

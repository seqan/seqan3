#include <cassert>
#include <string>

#include <seqan3/utility/type_traits/function_traits.hpp>

std::function my_caller = [] (size_t position, std::string & sequence)
{
    assert(position < sequence.size());
    return sequence[position];
};

using my_function_t = decltype(my_caller);

static_assert(std::same_as<seqan3::function_traits<my_function_t>::result_type, char>);
static_assert(seqan3::function_traits<my_function_t>::argument_count == 2);
static_assert(std::same_as<seqan3::function_traits<my_function_t>::argument_type_at<0>, size_t>);
static_assert(std::same_as<seqan3::function_traits<my_function_t>::argument_type_at<1>, std::string &>);

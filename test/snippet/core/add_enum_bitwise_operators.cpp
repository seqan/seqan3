#include <seqan3/core/add_enum_bitwise_operators.hpp>

using namespace seqan3;

int main() {}
//! [usage]
enum class my_enum
{
    VAL1 = 1,
    VAL2 = 2,
    COMB = 3
};

template <>
constexpr bool seqan3::add_enum_bitwise_operators<my_enum> = true;

my_enum e = my_enum::VAL1;
my_enum e2 = e | my_enum::VAL2;

// e2 == my_enum::COMB;
//! [usage]

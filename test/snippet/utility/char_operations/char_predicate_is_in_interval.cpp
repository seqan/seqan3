#include <seqan3/utility/char_operations/predicate.hpp>

int main()
{
    seqan3::is_in_interval<'A', 'G'>('C'); // returns true

    auto constexpr my_check = seqan3::is_in_interval<'A', 'G'>;
    my_check('H');  // returns false
}

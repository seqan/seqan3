#include <seqan3/utility/char_operations/predicate.hpp>

int main()
{
    seqan3::is_char<'C'>('C'); // returns true

    constexpr auto my_check = seqan3::is_char<'C'>;
    my_check('c'); // returns false, because case is different
}

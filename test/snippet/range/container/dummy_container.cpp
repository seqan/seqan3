#include <seqan3/range/container/dummy_container.hpp>
#include <seqan3/io/stream/debug_stream.hpp>

using namespace seqan3;

int main()
{
    dummy_container<char> dummy{"HALLO"}; // dummy has size 5 but stores no data

    debug_stream << dummy.size() << std::endl; // prints 5

    // debug_stream << dumm1[1] << std::endl; will NOT WORK, because there is no data to access

    // you can still modify the sequence, but only the size changes.
    dummy.insert(dummy.begin(), 5, 'O');  // now size == 10
    dummy.erase(dummy.begin());           // now size == 9

    debug_stream << dummy.size() << std::endl; // prints 9

    // as you saw above dummy.begin() gives you an iterator but you cannot dereference it:
    // debug_stream << *dumm1.begin() << std::endl; will NOT WORK, because there is no data to access
}

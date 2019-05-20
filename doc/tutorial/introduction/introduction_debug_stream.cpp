//! [debug]
#include <iostream>                          // for std::cerr, std::endl
#include <vector>                            // for std::vector
#include <seqan3/core/debug_stream.hpp> // for debug_stream

using namespace seqan3;

int main()
{
    std::vector<int> vec{-1,0,1};
    debug_stream << vec << std::endl;   // => [-1,0,1]
    // std::cerr << vec << std::endl;   // compiler error: no operator<< for std::vector<int>
}
//! [debug]

//! [debug]
#include <iostream>                          // for std::cerr, std::endl
#include <seqan3/core/debug_stream.hpp> 	// for debug_stream, includes <vector> as well

int main()
{
    std::vector<int> vec{-1,0,1};
    seqan3::debug_stream << vec << std::endl;   // => [-1,0,1]
    // std::cerr << vec << std::endl;   // compiler error: no operator<< for std::vector<int>
    return 0;
}
//! [debug]

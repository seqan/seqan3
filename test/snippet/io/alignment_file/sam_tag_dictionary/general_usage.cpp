//! [all]
#include <iostream>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/alignment_file/sam_tag_dictionary.hpp>

using namespace seqan3;

int main()
{
    sam_tag_dictionary dict{};         // initialise empty dictionary

    dict.get<"NM"_tag>() = 3;          // set SAM tag 'NM' to 3 (integer type)
    dict.get<"CO"_tag>() = "comment";  // set SAM tag 'CO' to "comment" (string type)

    auto nm = dict.get<"NM"_tag>();    // get SAM tag 'NM' (note: type is int32_t)
    auto co = dict.get<"CO"_tag>();    // get SAM tag 'CO' (note: type is std::string)

    debug_stream << nm << std::endl;      // will print '3'
    debug_stream << co << std::endl;      // will print "comment"
}
//! [all]

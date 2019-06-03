#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/view/convert.hpp>
#include <seqan3/range/view/to_rank.hpp>

using namespace seqan3;

int main()
{
//! [usage]
// This does not work:
// std::cout << 'C'_dna5;
// because the alphabet needs to be converted to char explicitly:
debug_stream << to_char('C'_dna5);  // prints 'C'

// The debug_stream, on the other hand, does this automatically:
debug_stream << 'C'_dna5;        // prints 'C'

// Vectors are also not printable to debug_stream:
std::vector<dna5> vec{"ACGT"_dna5};
//debug_stream << vec;
// but all types that model std::ranges::InputRange are printable to the debug_stream:
debug_stream << vec;            // prints "ACGT"

// ranges of non-alphabets are printed comma-separated:
debug_stream << (vec | view::to_rank); // prints "[0,1,2,3]"
//! [usage]

{
//! [flags]
uint8_t i = 71;
debug_stream << '\'' << i << "'\n";                                     // prints '71' (because flag is set by default)
debug_stream.unsetf(fmtflags2::small_int_as_number);                    // unsets the flag
debug_stream << '\'' << i << "'\n";                                     // prints 'G'
debug_stream << fmtflags2::small_int_as_number << '\'' << i << "'\n";   // prints '71' again
// instead of formatting the stream "inline", one can also call .setf()
//! [flags]
}

//! [set_underlying_stream]
{
    std::ostringstream o;
    debug_stream.set_underlying_stream(o);

    debug_stream << "ACGT"_dna5;

    o.flush();
    debug_stream << o.str(); // prints the string stream's buffer: "ACGT"
}

// this is UNDEFINED BEHAVIOUR, because the underlying stream went out-of-scope:
//debug_stream << 'C'_dna4;
//! [set_underlying_stream]

//! [set_underlying_stream2]
{
    std::ostringstream o;
    debug_stream_type my_stream{o};

    my_stream << "ACGT"_dna5;

    o.flush();
    debug_stream << o.str(); // prints the string stream's buffer: "ACGT"
}
// now your custom debug stream went out of scope with its underlying stream
//! [set_underlying_stream2]

}

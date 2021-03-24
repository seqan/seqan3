#if SEQAN3_HAS_ZLIB
//![example]
#include <seqan3/io/all.hpp>

// This does not work, the value must be overwritten within a function.
// seqan3::contrib::bgzf_thread_count = 1u;

int main()
{
    // Use one thread for (de-)compression.
    seqan3::contrib::bgzf_thread_count = 1u;

    // Read/Write compressed files.
    // ...
    return 0;
}
//![example]
#endif // SEQAN3_HAS_ZLIB

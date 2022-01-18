#if defined(SEQAN3_HAS_ZLIB)
//![example]
#include <seqan3/io/all.hpp>


// The `bgzf_thread_count` is a variable that can only be changed during run time of the program.
// This does not work, the value must be overwritten within a function.
// seqan3::contrib::bgzf_thread_count = 1u; // Doesn't work

int main()
{
    // Here we change the number of threads to `1`. This is a global
    // change and will effect the every future bgzf decompression call
    // which doesn't specify an explicit thread count parameter.
    // Already running decompression calls are unaffected by this
    seqan3::contrib::bgzf_thread_count = 1u;

    // Read/Write compressed files.
    // ...
    return 0;
}
//![example]
#endif // defined(SEQAN3_HAS_ZLIB)

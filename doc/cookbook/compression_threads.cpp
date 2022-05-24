#if defined(SEQAN3_HAS_ZLIB)
//![example]
#    include <seqan3/io/all.hpp>

// The `bgzf_thread_count` is a variable that can only be changed during the runtime of a program.
// The following does not work, the value must be overwritten within a function:
// seqan3::contrib::bgzf_thread_count = 1u; // Does not work.

int main()
{
    // Here, we change the number of threads to `1`.
    // This is a global change and will affect every future bgzf (de-)compression.
    // However, running (de-)compressions will not be affected.
    // `bgzf_thread_count` may be overwritten multiple times during the runtime of a program, in which case
    // the latest modification will determine the value.
    seqan3::contrib::bgzf_thread_count = 1u;

    // Read/Write compressed files.
    // ...
    return 0;
}
//![example]
#endif // defined(SEQAN3_HAS_ZLIB)

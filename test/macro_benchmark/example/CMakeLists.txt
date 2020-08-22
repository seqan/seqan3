
seqan3_macro_benchmark(
    memcpy_macro_benchmark.cpp
    REPETITIONS 10 # optional, default: 3
    META_DATA # optional
        foo bar # optional
        # seqan3_version 3.0.1 # automatically set
        # seqan3_git_version a1829272d7bafea48580af6eaac2bfbb678c3861 # automatically set
        # range-v3_version 0.10.0 # automatically set
        # range-v3_git_version 5daa38e80238a1b50a18c08b4a582178f414404c # automatically set
        # sdsl_version 3.0.0 # automatically set
        # sdsl_git_version 9a0d5676fd09fb8b52af214eca2d5809c9a32dbe # automatically set
)

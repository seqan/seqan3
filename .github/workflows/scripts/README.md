<!--
    SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
    SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
    SPDX-License-Identifier: CC-BY-4.0
-->

## g++_ram.sh

`VERSION` needs to be replaced with the gcc version. When `DO_TIME` is `0`, this is just a wrapper for gcc that does
nothing additionial. This mode might be necessary for CMake to properly detect ZLIB, etc..
When `DO_TIME` is `1`, every call will generate a `ram_usage.*` file with an unique extension containing the output
of `time -v`.

## process_compiler_error_log.py

Processes the output of compiling tests (`make`) and extracts errors.
Used in the API Stability and Latest Libraries workflow to create issues containing error messages.

## ram_usage.py

Processes a logfile containing all `time -v` outputs of `g++_ram.sh` and creates a table.
Used in the RAM Usage workflow.

An example input:

```
        Command being timed: "/usr/bin/g++-11 [...] -c /home/.../seqan3/test/unit/argument_parser/format_parse_test.cpp"
        User time (seconds): 9.21
        System time (seconds): 0.39
        Percent of CPU this job got: 99%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 0:09.64
        Average shared text size (kbytes): 0
        Average unshared data size (kbytes): 0
        Average stack size (kbytes): 0
        Average total size (kbytes): 0
        Maximum resident set size (kbytes): 672856
        Average resident set size (kbytes): 0
        Major (requiring I/O) page faults: 0
        Minor (reclaiming a frame) page faults: 211593
        Voluntary context switches: 9
        Involuntary context switches: 13
        Swaps: 0
        File system inputs: 8
        File system outputs: 44712
        Socket messages sent: 0
        Socket messages received: 0
        Signals delivered: 0
        Page size (bytes): 4096
        Exit status: 0
        Command being timed: "/usr/bin/g++-11 CMakeFiles/format_parse_test.dir/format_parse_test.cpp.o -o format_parse_test [...]"
        User time (seconds): 0.22
        System time (seconds): 0.04
        Percent of CPU this job got: 76%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 0:00.34
        Average shared text size (kbytes): 0
        Average unshared data size (kbytes): 0
        Average stack size (kbytes): 0
        Average total size (kbytes): 0
        Maximum resident set size (kbytes): 54872
        Average resident set size (kbytes): 0
        Major (requiring I/O) page faults: 0
        Minor (reclaiming a frame) page faults: 15478
        Voluntary context switches: 725
        Involuntary context switches: 0
        Swaps: 0
        File system inputs: 0
        File system outputs: 0
        Socket messages sent: 0
        Socket messages received: 0
        Signals delivered: 0
        Page size (bytes): 4096
        Exit status: 0
```

"Command being timed:" occurs every 23 lines.
The compiler is called two times: compiling and linking. Only the measurement for compiling is of interest.
The compile step will contain the path to the source file, e.g, `seqan3/test/unit/argument_parser/format_parse_test.cpp`.
The linking step will not, e.g., `CMakeFiles/format_parse_test.dir/format_parse_test.cpp.o -o format_parse_test`.

Hence, we iterate over all the time outputs (blocks of length 23) and check if there is `unit` in the command.
If yes, this is something we want the RAM usage of (`parsing_ram_usage`). If we want to parse the RAM usage, and we
are on the 9th line after the beginning of a block, we are in the line that contains `Maximum resident set size`.

The output is sorted by the RAM Usage (descending):
```
File,RAM in MiB
unit/search/search_test.cpp,2281
unit/alignment/pairwise/global_affine_unbanded_callback_test.cpp,1982
unit/search/search_collection_test.cpp,1620
unit/alignment/pairwise/align_pairwise_test.cpp,1492
unit/alignment/pairwise/alignment_configurator_test.cpp,1395
unit/alignment/pairwise/edit_distance/semi_global_edit_distance_max_errors_unbanded_test.cpp,1363
...
```

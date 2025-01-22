<!-- SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
     SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
     SPDX-License-Identifier: CC-BY-4.0
-->

SeqAn's pairwise sequence alignment algorithm is internally accelerated using multi-threading. The parallel execution
can be selected by specifying the seqan3::align_cfg::parallel configuration element. This will enable the asynchronous
execution of the alignments in the backend. For the user interface nothing changes as the returned
seqan3::algorithm_result_generator_range will preserve the order of the computed alignment results, i.e. the first
result corresponds to the first alignment given by the input range. By default, a thread pool with
std::thread::hardware_concurrency many threads will be created on a call to seqan3::align_pairwise and destructed when
all alignments have been processed and the seqan3::algorithm_result_generator_range goes out of scope. The configuration
element seqan3::align_cfg::parallel can be initialised with a custom thread count which determines the number of threads
that will be spawned in the background.<br>
Note that only independent alignment computations can be executed in parallel, i.e. you use this method when computing a
batch of alignments rather than executing them separately. <br>
Depending on your processor architecture you can gain a significant speed-up.

<!-- SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
     SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
     SPDX-License-Identifier: CC-BY-4.0
-->

 **Local Alignment** (better suited to find conserved segments):
```
                  tccCAGTTATGTCAGgggacacgagcatgcagagac
                     ||||||||||||
aattgccgccgtcgttttcagCAGTTATGTCAGatc
```
A \ref seqan3::align_cfg::method_local "local" alignment is effectively a global alignment of two partial sequences.
For example when two genes from different species are similar in short conserved regions and dissimilar in the
remaining regions. A global alignment would not find the local matching because it would try to align the entire
sequence. This is solved by the
[Smith-Waterman algorithm](https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm).

<!-- SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
     SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
     SPDX-License-Identifier: CC-BY-4.0
-->

| OP | Description                                                                                     |
|:--:|-------------------------------------------------------------------------------------------------|
| M  | Alignment match (can be a sequence match or mismatch, used only in basic CIGAR representations) |
| I  | Insertion to the reference                                                                      |
| D  | Deletion from the reference                                                                     |
| N  | Skipped region from the reference                                                               |
| S  | Soft clipping (clipped sequences present in seqan3::sam_record::sequence)                       |
| H  | Hard clipping (clipped sequences NOT present in seqan3::sam_record::sequence)                   |
| P  | Padding (silent deletion from padded reference)                                                 |
| =  | Sequence match                                                                                  |
| X  | Sequence mismatch

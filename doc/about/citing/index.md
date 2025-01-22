# Citing {#about_citing}

<!-- SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
     SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
     SPDX-License-Identifier: CC-BY-4.0
-->

[TOC]

SeqAn is a research project and depends strongly on being correctly attributed.
Please always cite the correct research papers when you develop software based
on SeqAn. This will help us continue to acquire funding and improve the library.

# Publications

SeqAn3 has not yet been published academically, for now cite SeqAn2:

| Publication   |                                                                                                   |
|---------------|---------------------------------------------------------------------------------------------------|
| Title         | The SeqAn C++ template library for efficient sequence analysis: A resource for programmers.       |
| Year          | 2017                                                                                              |
| Authors       | Reinert, K., Dadi, T. H., Ehrhardt, M., Hauswedell, H., Mehringer, S., Rahn, R., ... & Urgese, G. |
| Journal       | Journal of biotechnology, 261, 157-168                                                            |
| Links         | [original](https://doi.org/10.1016/j.jbiotec.2017.07.017)                                         |

Certain compontents of SeqAn are published separately. If you make strong use of
one of those compononts and/or specifically compare to that component, please
cite the respective publication **additionally**.

## Alignment module

| Publication   |                                                                                                   |
|---------------|---------------------------------------------------------------------------------------------------|
| Title         | Generic accelerated sequence alignment in SeqAn using vectorization and multi-threading.          |
| Year          | 2018                                                                                              |
| Authors       | Rahn, R., Budach, S., Costanza, P., Ehrhardt, M., Hancox, J., & Reinert, K                        |
| Journal       | Bioinformatics, 34(20), 3437-3445                                                                 |
| Links         | [original](https://doi.org/10.1093/bioinformatics/bty380)                                         |

Cite the above publication when you make strong use of the alignment module, in
particular if you rely on the high-performance computing capabilities.

## Search module

| Publication   |                                                                                                   |
|---------------|---------------------------------------------------------------------------------------------------|
| Title         | From theory to practice: Plug and play with succinct data structures.                             |
| Year          | 2014                                                                                              |
| Authors       | Gog, S., Beller, T., Moffat, A., & Petri, M.                                                      |
| Journal       | International Symposium on Experimental Algorithms (pp. 326-337). Springer, Cham.                 |
| Links         | [original](https://doi.org/10.1007/978-3-319-07959-2_28)                                          |


Full text indexing in SeqAn3 makes use of the Succint data structure library (SDSL).
Version 3 has not yet been published, for now cite SDSL-v2 (above).

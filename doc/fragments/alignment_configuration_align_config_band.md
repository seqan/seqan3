<!-- SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
     SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
     SPDX-License-Identifier: CC-BY-4.0
-->

SeqAn offers the computation of banded alignments to reduce the running time of the algorithm. This can be helpful if
the region in which the optimal alignment exists is known a priori. To specify the banded alignment the developer can
use the seqan3::align_cfg::band_fixed_size option.<br>
This band configuration is initialised with a seqan3::align_cfg::lower_diagonal and a seqan3::align_cfg::upper_diagonal.
The term diagonal is used to describe the position of the band boundary within the alignment matrix. The given value
represents the offset that the lower, respectively upper, diagonal is shifted from the main diagonal, which starts in
the origin of the alignment matrix. Accordingly, a negative value shifts the band boundary downwards in the alignment
matrix and a positive value shifts the band boundary to the right.
<br><br>
The band parameters might be restricted depending on the configured alignment algorithm, e.g. the origin of the
alignment matrix and the sink (the last cell in the last column) must be covered by the band when a global alignment is
ought to be computed.<br>
In general, the upper diagonal must always be greater than or equal to the lower diagonal to specify a valid band.

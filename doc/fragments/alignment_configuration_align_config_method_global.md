**Global Alignment**:
```
--T--CC-C-AGT--TATGT-CAGGGGACACG-A-GCATGCAGA-GAC
  |  || |  ||  | | | |||    || | | |  | ||||   |
AATTGCCGCC-GTCGT-T-TTCAG----CA-GTTATG-T-CAGAT--C
```
Finding the optimal global alignment of two sequences is solved by the
[Needleman-Wunsch algorithm](https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm).<br>

**Semi-global Alignment** (e.g. overlapping sequences):
```
                  TCCCAGTTATGTCAGgggacacgagcatgcagagac
                  |||||||||||||||
aattgccgccgtcgttttTCCCAGTTATGTCAG
```
The semi-global alignment is a specially configured global alignment, namely we do not penalize gaps at the ends of
the alignment. Semi-global alignments are often used in genome assembly applications when trying to find matching
overlaps.

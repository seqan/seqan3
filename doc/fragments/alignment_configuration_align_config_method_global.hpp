namespace seqan3::doxygen
{
/*!\details
 *
 * **Global Alignment**:\verbatim
--T--CC-C-AGT--TATGT-CAGGGGACACG-A-GCATGCAGA-GAC
  |  || |  ||  | | | |||    || | | |  | ||||   |
AATTGCCGCC-GTCGT-T-TTCAG----CA-GTTATG-T-CAGAT--C\endverbatim
 * Finding the optimal global alignment of two sequences is solved by the
 * [Needleman-Wunsch algorithm](https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm).
 *
 * **Semi-global Alignment** (e.g. overlapping sequences):
 * \verbatim
                  TCCCAGTTATGTCAGgggacacgagcatgcagagac
                  |||||||||||||||
aattgccgccgtcgttttTCCCAGTTATGTCAG
\endverbatim
 * The semi-global alignment is a specially configured global alignment, namely we do not penalize gaps at the ends of
 * the alignment. Semi-global alignments are often used in genome assembly applications when trying to find matching
 * overlaps.
 */
using alignment_configuration_align_config_method_global = void;
}

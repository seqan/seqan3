namespace seqan3::doxygen
{
/*!\details
 *
 * **Local Alignment** (better suited to find conserved segments):
 * \verbatim
                  tccCAGTTATGTCAGgggacacgagcatgcagagac
                     ||||||||||||
aattgccgccgtcgttttcagCAGTTATGTCAGatc
\endverbatim
 * A \ref seqan3::align_cfg::method_local "local" alignment is effectively a global alignment of two partial sequences.
 * For example when two genes from different species are similar in short conserved regions and dissimilar in the
 * remaining regions. A global alignment would not find the local matching because it would try to align the entire
 * sequence. This is solved by the
 * [Smith-Waterman algorithm](https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm).
 */
using alignment_configuration_align_config_method_local = void;
}

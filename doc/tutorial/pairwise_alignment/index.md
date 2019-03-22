# Pairwise Alignment {#tutorial_pairwise_alignment}

<b>Learning Objective:</b> <br/>

In this tutorial you will learn how to compute pairwise sequence alignments with SeqAn3.
This tutorial is a walkthrough with links to the API documentation and is also meant as a source for copy-and-paste code.

\tutorial_head{Intermediate, 60-90 min, \ref setup \ref tutorial_alphabets \ref tutorial_ranges, }

[TOC]

---

# Introduction

Aligning biological sequences is a very prominent component in many bioinformatics applications and pipelines.
Well known genomic applications that use pairwise sequence alignments are read mapping, genome assembly, variant
detection, multiple sequence alignment as well as protein search.

The goal of the pairwise alignment is to obtain an optimal transcript that describes how two sequences are related to
each other by means of substitutions, insertions, or deletions. The computed transcript describes then the operations
necessary to translate the one sequence into the other, as can be seen in the following picture.

\htmlonly
<img src="doc/tutorial/pairwise_alignment/align_transcript.png" alt="A transcript between two aligned sequences" style="width:300px;"/>
\endhtmlonly

The alignment problem is solved with a dynamic programming (DP) algorithm which runs in \f$ (\mathcal{O}(n^2))\f$ time
and space. Besides the global alignment approach many more variations of this DP based algorithm have been developed
over time. SeqAn unified all of these approaches into a single DP core implementation which can be extended easily and
thus, with all possible configurations, is a very versatile and powerful tool to compute many desired alignment variants.

# Computing pairwise alignments

Let us first have look at an example of computing a global alignment in SeqAn3.

\includelineno doc/tutorial/pairwise_alignment/simple_global_alignment.cpp

In the above example we want to compute the similarity of two seqan3::dna4 sequences using a global alignment.
For this we need to import the dna4 alphabet, the seqan3::nucleotide_scoring_scheme and the seqan3::align_pairwise in
the beginning of our file. We also import the seqan3::debug_stream which allows us to print various types in a nice
formatted manner.

In the beginning of the file we are defining our two sequences `first_seq` and `second_seq`.
If you feel puzzled about what the `_dna4` suffix does we recommend to read the \ref tutorial_alphabets and
\ref tutorial_ranges before which gives a good introduction to this topic.
In line 16-17 we configure the alignment job with the most simplistic configuration possible.
In this case it is a global alignment with edit distance.
Later in this tutorial we will give a more detailed description of the \ref alignment_configurations "configuration" and
how it can be used.
The minimum requirement for computing a pairwise alignment is to specify the seqan3::align_cfg::mode and the
seqan3::align_cfg::scoring configuration elements. The first one selects the internal algorithm and the second one
provides the scoring scheme that should be used to score a pair of sequence characters.

Now we are going to call seqan3::align_pairwise. This interface requires two arguments: a tuple or a range of tuples
with exactly two elements and the configuration object. Independent of the number of pairwise alignment jobs submitted
to the algorithm it always returns a range over seqan3::alignment_result. Later we will see how we can use a
continuation interface which calls a user-defined function rather than iterating over the results sequentially.
Finally we output the score for the computed alignment.

\attention Just calling `pairwise_align` does nothing as it returns a range which is evaluated in a lazy manner.
           Only when calling begin or incrementing the iterator over the range the alignment computation is invoked.

\assignment{Assignment 1}
Copy and paste the minimal working example from above into a cpp file in your working directory and compile and run it.
Add the two sequences "ACGTGACTGACT" and "AGGTACGAGCGACACT" to the set and compute all pairwise sequence alignments
of the four sequences and output the results.

\endassignment
\solution

\include doc/tutorial/pairwise_alignment/simple_global_alignment_solution_1.cpp

First we create the vector of seqan3::dna4 sequences. We keep the configuration as is and then modify the initial code
to a range-based for loop looping over the alignment results. Since the seqan3::alignment_result is a class template and the
template parameters are determined during the configuration step we use auto as the result type.
The current result is cached inside of the lazy range and we capture the result as `const &` in order to not tamper with
the result values.

\endsolution

Congratulations, you have computed your first pairwise alignment with seqan3!
As you can see, the interface is really simple, yet the configuration object makes it extremely flexible to conduct
various different alignment calculations. In the following chapter you will learn more about the various configurations.

# Alignment configurations
\anchor alignment_configurations

The configuration object is the core of the alignment interface. It allows to easily configure the alignment algorithm
without changing the interface of the pairwise alignment. It uses a seqan3::configuration object to chain different
configuration elements together using the logical or-operator ('|'-operator).
You can find an overview over the available configurations \ref configuration "here".

The configuration elements are all classes that wrap the actual information necessary for the configuration of the
alignment algorithm. Depending on the configuration specification certain features of the algorithm are enabled or
disabled. Moreover, during the initialisation of the algorithm the best implementation is chosen based on the given
configurations. To avoid possible ambiguities with the configurations of other algorithms, the configuration elements
for the alignment live in the special namespace seqan3::align_cfg.

## Global and semi-global alignment

The most straightforward algorithm is the global alignment which can be configured using the seqan3::align_cfg::mode
together with seqan3::global_alignment as the constructor argument to the configuration element.
The seqan3::align_cfg::mode must be given as a minimal configuration together with the seqan3::align_cfg::scoring.
The global alignment can be further refined by setting the seqan3::align_cfg::aligned_ends option.
The aligned ends specify wether or not gaps at the end of the sequences are penalised.
In SeqAn you can configure this behaviour for every end (front and back of the first sequence and second sequence)
separately using the seqan3::align_config::end_gaps class. However, there are \ref predefined_end_gap_configurations
"predefined" sets that allow you to swiftly choose the option that fits your needs best.

\assignment{Assignment 2}
Use the sequences from above and compute a semi-global alignment where the ends of the second sequence are free.

\endassignment
\solution

\include doc/tutorial/pairwise_alignment/pa_assignment_2_solution.cpp

\endsolution

## Scoring schemes

The seqan3::align_cfg::scoring configures the scoring scheme used to score two alphabet values. It must be constructed
with a class that models the seqan3::scoring_scheme_concept. In SeqAn we have two scoring schemes available: the
seqan3::nucleotide_scoring_scheme and the seqan3::amino_acid_scoring_scheme. Use the first one when you align
nucleotide sequences and the latter when working with amino acid sequences.

\assignment{Assignment 3}
Make the assignment from above work with protein sequences using BLOSUM62 matrix.

\endassignment
\solution

\include doc/tutorial/pairwise_alignment/pa_assignment_3_solution.cpp

\endsolution

## Gap schemes

The seqan3::align_cfg::gap configures the gap scheme. The gap scheme determines how gaps are penalised. If no
gap scheme is provided to the alignment configuration it will automatically choose -1 for a gap extension and 0 for a
gap opening.

\note In SeqAn3 the gap open score is an additional cost to the gap score that is added when a gap is opened in affine
alignment algorithm.

\todo assignment compute protein alingment with affine gaps

## Alignment result

The seqan3::align_cfg::result regulates the outcome of the alignment algorithm. In the default case only the score
will be computed. To obtain also the alignment or the begin and end coordinates of the alignment the
seqan3::align_cfg::result has to be configured accordingly.
Depending on the configuration the algorithm choses the most efficient implementation to obtain the requested result.
The result can be finally examined using the seqan3::alignment_result class.

\todo repeat the protein alignment and print the alignment usign the debug stream.

## Edit distance

There is also a shortcut for computing the edit distance for two sequences. The edit distance is a special metric
to count the number of edit operations between two sequences. It can be solved with a fast bit-vector algorithm.
To compute the edit-distance you can use the shortcut seqan3::align_cfg::edit. This only works for nucleotide
sequences. It can be further refined using seqan3::align_cfg::end_gaps::free_ends_first and the
seqan3::align_cfg::max_error configuration.

\todo Compute all pairwise edit operations of the sequences and filter out every alignment with a score higher than 5

## Invalid configurations

Since the configuration allows to possibly chain any configuration element together we added some security details.
First, not all configuration elements of the alignment can be combined with each other. In this case a static assert
will inform you about the invalid combination of configurations. You can find \ref configuration "here " a table with
the information about the valid configurations.

<!-- # Parallel execution -->

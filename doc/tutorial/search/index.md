# Indexing and searching with SeqAn {#tutorial_index_search}

<b>Learning Objective:</b><br>
In this tutorial you will learn how to construct an index and conduct searches.

\tutorial_head{Moderate,
              60 Minutes,
              \ref tutorial_ranges (Recommended)<br>\ref tutorial_pairwise_alignment (only for last assignment),
              [FM-Index paper](https://doi.org/10.1109/SFCS.2000.892127)<br>
              [FM-Index on Wikipedia](https://en.wikipedia.org/wiki/FM-index)<br>
              [Bi-FM-Index paper](https://doi.org/10.1109/BIBM.2009.42)}

[TOC]

---

# Introduction

Exact and approximate string matching is a common problem in bioinformatics, e.g. in read mapping.
Usually, we want to search for a match of one or many small sequences (queries) in a database consisting of one or
more large sequences (references). Trivial approaches such as on-line search fail at this task for large data due
to prohibitive run times.

The general solution is to use a data structure called **index**. An index is built over the reference and only needs
to be computed once for a given data set. It is used to speed up the search in the reference and can be re-used by
storing it to disk and loading it again without recomputation.

There are two groups of indices: hash tables and suffix-based indices. SeqAn implements the FM-Index and Bi-FM-Index
as suffix-based indices. While the Bi-FM-Index needs more space (a factor of about 2), it allows faster searches.

Given an index, SeqAn will choose the best search strategy for you. Since very different algorithms may be selected
internally depending on the configuration, it is advisable to do benchmarks with your application. A rule of thumb is
to use the Bi-FM-Index when allowing more than 2 errors.

This tutorial will show you how to use the seqan3::fm_index and seqan3::bi_fm_index to create indices and how to
search them efficiently using seqan3::search.

## Capabilities

With this module you can:
* Create, store and load (Bi)-FM-Indices
* Search for exact matches
* Search for approximate matches (allowing substitutions and indels)

The results of the search can be passed on to other modules, e.g. to create an alignment.

## Terminology
### Reference
A reference is the data you want to search in, e.g. a genome or protein database.
### Query
A query is the data you want to search for in the reference, e.g. a read or an amino acid sequence.
### Index
An index is a data structure built over the reference that allows fast searches.
### FM-Index
The full-text minute index (FM-Index) is an index that is similar to a suffix tree, but much smaller.
It is used by most state-of-the-art read mappers and aligners. You can find more information on FM-Indicies
[in the original publication](https://doi.org/10.1109/SFCS.2000.892127) and
[on Wikipedia](https://en.wikipedia.org/wiki/FM-index).
### Bi-FM-Index
The bidirectional FM-Index (Bi-FM-Index) is an extension of the FM-Index that enables faster searches, especially
when allowing multiple errors. But it uses almost twice the amount of memory the FM-Index uses. You can find
more information on Bi-FM-Indicies [here](https://doi.org/10.1109/BIBM.2009.42).

## Example

Constructing a (Bi-)FM-Index is very simple:

\include doc/tutorial/search/search_basic_index.cpp

You can also index text collections (e.g. genomes with multiple chromosomes or protein databases):

\snippet doc/tutorial/search/search_small_snippets.cpp text_collection

The indices can also be stored and loaded from disk by using cereal.

\snippet doc/tutorial/search/search_small_snippets.cpp store

\snippet doc/tutorial/search/search_small_snippets.cpp load

Note that in contrast to the construction via a given `text`, the template cannot be deduced by the compiler when
using the default constructor so you have to provide template arguments.
<a name="assignment_create_index"></a>
\assignment{Assignment 1}
You are given the text
\code
dna4_vector text{"CGCTGTCTGAAGGATGAGTGTCAGCCAGTGTAACCCGATGAGCTACCCAGTAGTCGAACTGGGCCAGACAACCCGGCGCTAATGCACTCA"_dna4};
\endcode
Create a seqan3::fm_index over the reference, store the index and load the index into a new seqan3::fm_index object.
Print whether the indices are identical or differ.
\endassignment

\solution
\snippet doc/tutorial/search/search_solution1.cpp solution
**Expected output:**
```console
The indices are identical!
```
\endsolution

# Search

Using an index, we can now conduct searches for a given query. In this part we will learn how to search exactly,
allow substitutions and indels, and how to configure what kind of results we want, e.g. all results vs. only
the best result.

## Terminology
### Match
If a query can be found in the reference with the specified error configuration, we denote this result as a match.
If the query matches the reference without any errors, we call this an exact match. If a match contains errors,
it's an approximate match.

## Searching for exact matches

We can search for all exact matches using seqan3::search:

\include doc/tutorial/search/search_basic_search.cpp

When using a single text, the search will return a std::vector of the begin positions of matches in the reference.

For text collections, a std::vector over std::tuple where the first element represents the index of the text in the
collection and the second element represents the position within the text is returned.

You can also pass multiple queries at the same time:

\snippet doc/tutorial/search/search_small_snippets.cpp multiple_queries

The returned result will be a vector of individual results. Each element refers to the respective query.
<a name="assignment_exact_search"></a>
\assignment{Assignment 2}
Search for all exact occurrences of `GCT` in the text from [assignment 1](#assignment_create_index).<br>
Print the number of hits and the positions in a ascending ordering.<br>
Do the same for the following text collection:
\code
std::vector<dna4_vector> text{"CGCTGTCTGAAGGATGAGTGTCAGCCAGTGTA"_dna4,
                              "ACCCGATGAGCTACCCAGTAGTCGAACTG"_dna4,
                              "GGCCAGACAACCCGGCGCTAATGCACTCA"_dna4};
\endcode
\endassignment

\solution
\include doc/tutorial/search/search_solution2.cpp
**Expected output:**
```console
=====   Running on a single text   =====
There are 3 hits.
The positions are [1,41,77]

===== Running on a text collection =====
There are 3 hits.
The positions are [(0,1),(1,9),(2,16)]
```
\endsolution

## Searching for approximate matches

If you want to allow for errors in your query, you need to configure the approximate search with a
seqan3::search_cfg::max_error or seqan3::search_cfg::max_error_rate object.
These are constructed with four parameters:

With absolute numbers (seqan3::search_cfg::max_error):

- seqan3::search_cfg::total: Maximum number of total errors
- seqan3::search_cfg::substitution: Maximum number of substitutions
- seqan3::search_cfg::insertion: Maximum number of insertions
- seqan3::search_cfg::deletion: Maximum number of deletions

Use seqan3::search_cfg::max_error_rate with the same elements to define error rates \f$\in[0,1]\f$.

\note When using >= 2 errors it is advisable to use a Bi-FM-Index since searches will be faster.

To search for either 1 insertion or 1 deletion you can write:

\snippet doc/tutorial/search/search_small_snippets.cpp error_search

Note that the error amount/rates do not have to sum up to the total of errors allowed:
\snippet doc/tutorial/search/search_small_snippets.cpp error_sum

Allows 2 errors which can be any combination of 2 substituions, 1 insertion and 1 deletion.

Defining only the total will set all error types to this value.

Defining any of the error types, but no total, will set the total to the sum of all error types; if a total is set,
it will act as an upper bound.
<a name="assignment_approximate_search"></a>
\assignment{Assignment 3}
Search for all occurrences of `GCT` in the text from [assignment 1](#assignment_create_index).

Allow up to 1 substitution and print the matching substring of the reference.
\hint The search will give you positions in the text. To access the corresponding subrange of the text you can use
std::span:
\include doc/tutorial/search/search_span.cpp
\endhint
\endassignment

\solution
\include doc/tutorial/search/search_solution3.cpp
**Expected output:**
```console
There are 12 hits.
At position 1: GCT
At position 5: TCT
At position 12: GAT
At position 23: GCC
At position 36: GAT
At position 41: GCT
At position 57: ACT
At position 62: GCC
At position 75: GCG
At position 77: GCT
At position 83: GCA
At position 85: ACT
```
\endsolution

## Search modes

Besides the error configuration, you can define what kind of hits should be reported:

- seqan3::search_cfg::all: Report all hits that satisfy the (approximate) search.
- seqan3::search_cfg::best: Report the best hit, i.e. the *first* hit with the lowest edit distance.
- seqan3::search_cfg::all_best: Report all hits with the lowest edit distance.
- seqan3::search_cfg::strata: best+x mode. Report all hits within the x-neighbourhood of the best hit.

The mode is appended to the error configuration by using the `|`-operator:
\snippet doc/tutorial/search/search_small_snippets.cpp mode_best

The strata mode needs an additional parameter:
\snippet doc/tutorial/search/search_small_snippets.cpp mode_strata

If the best hit had an edit distance of 1, the strata mode would report all hits with up to an edit distance of 3.
Since in this example the total error number is set to 2, all hits with 1 or 2 errors would be reported.

\assignment{Assignment 4}
Search for all occurrences of `GCT` in the text from [assignment 1](#assignment_create_index).<br>
Allow up to 1 error of any type and print the number of hits for each search mode (use `mode::strata{1}`).
\endassignment

\solution
\include doc/tutorial/search/search_solution4.cpp
**Expected output:**
```console
Searching all hits
There are 25 hits.
Searching all best hits
There are 3 hits.
Searching best hit
There is 1 hit.
Searching all hits in the 1-stratum
There are 25 hits.
```
\endsolution

## One last exercise

\assignment{Assignment 5}
Search for all occurrences of `GCT` in the text from [assignment 1](#assignment_create_index).<br>
Allow up to 1 error of any type search for all occurrences in the all_best mode.<br>
Align the query to each of the found positions in the genome and print the score and alignment.<br>
**BONUS**<br>
Do the same for the text collection from [assignment 2](#assignment_exact_search).
\hint Similarly to [assignment 3](#assignment_approximate_search), you have to extract the subrange of the text in
order to align against it:
```cpp
// std::span returnes the subrange [start, start + span - 1], i.e. span many elements beginning at start

// For single texts.
std::span text_view{std::data(text) + start, span};

// For text collections. idx is the current index of the text collection returned by the search.
std::span text_view{std::data(text[idx]) + start, span};
```
\endhint
\endassignment

\solution
\include doc/tutorial/search/search_solution5.cpp
**Expected output:**
```console
Searching all best hits allowing for 1 error in a single text
There are 3 hits.
-----------------
score:    0
database: GCT
query:    GCT
=============
score:    0
database: GCT
query:    GCT
=============
score:    0
database: GCT
query:    GCT
=============

Searching all best hits allowing for 1 error in a text collection
There are 3 hits.
-----------------
score:    0
database: GCT
query:    GCT
=============
score:    0
database: GCT
query:    GCT
=============
score:    0
database: GCT
query:    GCT
=============
```
\endsolution

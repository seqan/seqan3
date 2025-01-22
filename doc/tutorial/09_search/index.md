# Indexing and searching with SeqAn {#tutorial_index_search}

<!-- SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
     SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
     SPDX-License-Identifier: CC-BY-4.0
-->

<b>Learning Objective:</b><br>
In this tutorial, you will learn how to construct an index and conduct searches.

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
Usually, we want to search for a hit of one or many small sequences (queries) in a database consisting of one or
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
* Search for exact hits
* Search for approximate hits (allowing substitutions and indels)

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

\include doc/tutorial/09_search/search_basic_index.cpp

You can also index text collections (e.g. genomes with multiple chromosomes or protein databases):

\snippet doc/tutorial/09_search/search_small_snippets.cpp text_collection

The indices can also be stored and loaded from disk by using cereal.

\snippet doc/tutorial/09_search/search_small_snippets.cpp store

\snippet doc/tutorial/09_search/search_small_snippets.cpp load

Note that in contrast to the construction via a given `text`, the template cannot be deduced by the compiler when
using the default constructor so you have to provide template arguments.
\anchor assignment_create_index
\assignment{Assignment 1}
You are given the text
\code
dna4_vector text{"CGCTGTCTGAAGGATGAGTGTCAGCCAGTGTAACCCGATGAGCTACCCAGTAGTCGAACTGGGCCAGACAACCCGGCGCTAATGCACTCA"_dna4};
\endcode
Create a seqan3::fm_index over the reference, store the index and load the index into a new seqan3::fm_index object.
Print whether the indices are identical or differ.
\endassignment

\solution
\snippet doc/tutorial/09_search/search_solution1.cpp solution
**Expected output:**
```console
The indices are identical!
```
\endsolution

# Search

Using an index, we can now conduct searches for a given query. In this part, we will learn how to search exactly,
allow substitutions and indels, and how to configure what kind of results we want, e.g. all results vs. only
the best result.

## Terminology

**Exact search:** Finds all locations of the query in the reference without any errors.

**Approximate search:** Finds all locations of the query in the reference with substitutions and indels within the confined set of maximal allowed errors.

**Hit:** A single result that identifies a specific location in the reference where a particular query can be found by either an exact or an approximate search.

## Searching for exact hits

We can search for all exact hits using seqan3::search:

\include doc/tutorial/09_search/search_basic_search.cpp

You can also pass multiple queries at the same time:

\snippet doc/tutorial/09_search/search_small_snippets.cpp multiple_queries

The returned result is a lazy range over individual results, where each entry represents a specific location
within the reference sequence for a particular query.

\anchor assignment_exact_search
\assignment{Assignment 2}
Search for all exact occurrences of `GCT` in the text from [assignment 1](#assignment_create_index).<br>
Print the number of hits and their positions within the reference sequence.<br>
Do the same for the following text collection:
\code
std::vector<dna4_vector> text{"CGCTGTCTGAAGGATGAGTGTCAGCCAGTGTA"_dna4,
                              "ACCCGATGAGCTACCCAGTAGTCGAACTG"_dna4,
                              "GGCCAGACAACCCGGCGCTAATGCACTCA"_dna4};
\endcode
\endassignment

\solution
\include doc/tutorial/09_search/search_solution2.cpp
**Expected output:**
```console
=====   Running on a single text   =====
There are 3 hits.
The following hits were found:
<query_id:0, reference_id:0, reference_pos:1>
<query_id:0, reference_id:0, reference_pos:41>
<query_id:0, reference_id:0, reference_pos:77>

===== Running on a text collection =====
There are 3 hits.
The following hits were found:
<query_id:0, reference_id:0, reference_pos:1>
<query_id:0, reference_id:1, reference_pos:9>
<query_id:0, reference_id:2, reference_pos:16>
```
\endsolution

## Searching for approximate hits

Up until now, we have seen that we can call the search with the query sequences and the index. In addition, we
can provide a third parameter to provide a user defined
\ref search_configuration_section_introduction "search configuration".
If we do not provide a user defined search configuration, the seqan3::search_cfg::default_configuration will be used,
which triggers an exact search, finding all hits for a particular query. In the following, we will see how we can
change the behaviour of the search algorithm by providing a user defined search configuration.

### Max error configuration

You can specify the error configuration for the approximate search using the
\ref search_configuration_subsection_error "seqan3::search_cfg::max_error_*" configuration.
Here you can use a combination of the following configuration elements to specify exactly how the errors can be
distributed during the search:
 - seqan3::search_cfg::max_error_total,
 - seqan3::search_cfg::max_error_substitution,
 - seqan3::search_cfg::max_error_insertion and
 - seqan3::search_cfg::max_error_deletion.

Each of the configuration elements can be constructed with either an absolute number of errors or an error rate
depending on the context. These are represented by the following types:

- seqan3::search_cfg::error_count: Absolute number of errors
- seqan3::search_cfg::error_rate: Rate of errors \f$\in[0,1]\f$

By combining the different error types using the `|`-operator, we give you full control over the error distribution.
Thus, it is possible to set an upper limit of allowed errors but also to refine the error distribution by specifying the
allowed errors for substitutions, insertions, or deletions.

\note When using >= 2 errors it is advisable to use a Bi-FM-Index since searches will be faster.

For example, to search for either 1 insertion or 1 deletion you can write:

\snippet doc/tutorial/09_search/search_small_snippets.cpp error_search

Here, we restrict the approximate search to only allow one error. This can be then either an insertion or a deletion
but not both, since that would exceed the total error limit.
This basically means that the error counts/rates do not have to sum up to the total of errors allowed:

\snippet doc/tutorial/09_search/search_small_snippets.cpp error_sum

In the above example, we allow 2 errors, which can be any combination of 2 substitutions, 1 insertion and 1 deletion.
Defining only the total will set all error types to this value, i.e. if the total error is set to an error count of 2, any
combination of 2 substitutions, 2 insertions and 2 deletions is allowed.
On the other hand, when defining any of the error types but no total, the total will be set to the sum of all error
types. For example, if we would not specify a total error of 1 in the first example above, the total error would be set
to 2 automatically. Hence, the search will also include approximate hits containing one insertion and one deletion.

<a name="assignment_approximate_search"></a>
\assignment{Assignment 3}
Search for all occurrences of `GCT` in the text from [assignment 1](#assignment_create_index).

Allow up to 1 substitution and print all occurrences.
\endassignment

\solution
\include doc/tutorial/09_search/search_solution3.cpp
**Expected output:**
```console
[<query_id:0, reference_id:0, reference_pos:1>,
 <query_id:0, reference_id:0, reference_pos:5>,
 <query_id:0, reference_id:0, reference_pos:12>,
 <query_id:0, reference_id:0, reference_pos:23>,
 <query_id:0, reference_id:0, reference_pos:36>,
 <query_id:0, reference_id:0, reference_pos:41>,
 <query_id:0, reference_id:0, reference_pos:57>,
 <query_id:0, reference_id:0, reference_pos:62>,
 <query_id:0, reference_id:0, reference_pos:75>,
 <query_id:0, reference_id:0, reference_pos:77>,
 <query_id:0, reference_id:0, reference_pos:83>,
 <query_id:0, reference_id:0, reference_pos:85>]
```
\endsolution

## Which hits are reported?

Besides the max error configuration, you can specify the scope of the search algorithm. This means that you can control
which hits should be reported by the search.

### Hit configuration

To do so, you can use one of the following \ref search_configuration_subsection_hit_strategy "seqan3::search_cfg::hit_*"
configurations:

- seqan3::search_cfg::hit_all: Report all hits that satisfy the (approximate) search.
- seqan3::search_cfg::hit_single_best: Report the best hit, i.e. the *first* hit with the lowest edit distance.
- seqan3::search_cfg::hit_all_best: Report all hits with the lowest edit distance.
- seqan3::search_cfg::hit_strata: best+x strategy. Report all hits within the x-neighbourhood of the best hit.

In contrast to the max error configuration, which allows a combination of the different error configuration objects, the hit
configuration can only exist once within one search configuration. Trying to specify more than one hit configuration
in one search configuration will fail at compile time with a static assertion.
Sometimes the program you write requires to choose between different hit configurations depending on a user given
program argument at runtime. To handle such cases you can also use the dynamic configuration seqan3::search_cfg::hit.
This configuration object represents one of the four hit configurations mentioned previously and can be modified at
runtime. The following snippet gives an example for this scenario:

\snippet doc/tutorial/09_search/search_small_snippets.cpp hit_dynamic

Note that the same rules apply to both the dynamic and static hit configuration. That is, it can be added via the
`|`-operator to the search configuration but cannot be combined with any other hit configuration.

A closer look at the strata configuration reveals that it is initialised with an additional parameter called the
`stratum`. The `stratum` can be modified even after it was added to the search configuration like the following example
demonstrates:

\snippet doc/tutorial/09_search/search_small_snippets.cpp hit_strata

Here we introduced a new concept when working with the seqan3::configuration object, which is much like the access
interface of a std::tuple. Concretely, it is possible to access the stored
configuration using the `get<cfg_type>(cfg)` interface, where `cfg_type` is the name of the configuration type we would like
to access. The `get` interface returns a reference to the stored object that is identified by the given name. If you
try to access an object which does not exist within the search configuration, a static assert will be emitted at compile
time such that no invalid code can be generated.

\note We need to use the expression `using seqan3::get;` before we can call the `get` interface in order to allow the
      compiler to find the correct implementation of it based on the passed argument. This is related to how C++
      resolves unqualified lookup of free functions in combination with function templates using an explicit template
      argument such as the `get` interface does.

So, the open question remains what the stratum actually does. In the above example, if the best hit found by
the search for a particular query had an edit distance of 1, the strata strategy would report all hits with up to an
edit distance of 2.
Since in this example the total error number is set to 2, all hits with 1 or 2 errors would be reported.

\assignment{Assignment 4}
Search for all occurrences of `GCT` in the text from [assignment 1](#assignment_create_index).<br>
Allow up to 1 error of any type and print the number of hits for each hit strategy (use
`seqan3::search_cfg::strata{1}`).
\hint
    You can use std::ranges::distance to get the size of any range. Depending on the underlying range properties, this
    algorithm will use the optimal way to compute the number of elements contained in the range.
\endhint
\endassignment

\solution
\include doc/tutorial/09_search/search_solution4.cpp
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

## Controlling the search output

When calling the search algorithm, a lazy range over seqan3::search_result objects is returned. Each result object represents a
single hit. This means that merely calling the seqan3::search algorithm will do nothing except configure the search
algorithm based on the given search configuration, query and index. Only when iterating over the lazy search result
range, the actual search for every query is triggered. We have done this automatically in the previous examples when
printing the result to the seqan3::debug_stream which then invokes a range based iteration over the returned range or
by using the std::ranges::distance algorithm. However, in many cases, we want to access the specific positions and
information stored in the seqan3::search_result object to proceed with our application. Since some information might be
more compute-intensive than others, there is a way to control what the final search result object will contain.

### Output configuration

The behaviour of the search algorithm is further controlled through the
\ref search_configuration_subsection_output "seqan3::search_cfg::output_*" configurations.
The following output configurations exists:

- seqan3::search_cfg::output_query_id
- seqan3::search_cfg::output_reference_id
- seqan3::search_cfg::output_reference_begin_position
- seqan3::search_cfg::output_index_cursor

Similarly to the max error configurations, you can arbitrarily combine the configurations to customise the final output.
For example, if you are only interested in the position of the hit within the reference sequence, you can use the
seqan3::search_cfg::output_reference_begin_position configuration. Instead, if you need access to the index where the
hit was found, you can use the seqan3::search_cfg::output_index_cursor configuration.

\note If you do not provide any output configuration, then the query id and reference id as well as the
      reference begin position will be automatically reported. If you select only one in your search configuration, then only this one
      will be available in the final search result.

## One last exercise

In the final example, we will extend our previous search examples to also compute the alignment of the found hits and their
respective reference infixes. To do so, we recommend working through the \ref tutorial_pairwise_alignment tutorial
first.

\assignment{Assignment 5}
Search for all occurrences of `GCT` in the text from [assignment 1](#assignment_create_index).<br>
Allow up to 1 error of any type and search for all occurrences with the strategy `seqan3::search_cfg::hit_all_best`.<br>
Align the query to each of the found positions in the genome and print the score and alignment.<br>
**BONUS**<br>
Do the same for the text collection from [assignment 2](#assignment_exact_search).

\hint
The search will give you positions in the text. To access the corresponding subrange of the text you can use
std::span:
\include doc/tutorial/09_search/search_span.cpp
\endhint
\endassignment

\solution
\include doc/tutorial/09_search/search_solution5.cpp
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

# Minimisers {#tutorial_minimiser}

<!-- SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
     SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
     SPDX-License-Identifier: CC-BY-4.0
-->

[TOC]

This tutorial introduces minimisers. Minimisers are a compact representation of a DNA or RNA sequence that is closely
related to, but more efficient than, a representation by k-mers. Minimisers are implemented as a view in SeqAn3. For
more information about views and how to implement your own view, please see
\ref tutorial_ranges and \ref howto_write_a_view .

\tutorial_head{Easy, 20 min, , \ref tutorial_ranges\, \ref howto_write_a_view}

# Motivation

A common way to work with sequences is to obtain their k-mers, but even with a large size of k, the number of k-mers is
non-trivial, therefore storing k-mers is costly. At the same time, due to the overlap of consecutive k-mers, the
information they contain is highly redundant. That is why minimisers come in handy. Minimisers are k-mers, which have
a minimal value in a window of a specified size.  Minimal could for example mean lexicographically smallest. By storing
only these minimal k-mers the storage cost is significantly reduced while maintaining a similar amount of information.

# Minimiser Workflow

Because minimisers are minimal k-mers, they depend on a given k-mer size, a given shape (specifying which positions
should be considered) and a window size, which has to be greater than or equal to the k-mer size. If all these values are
given, then the minimisers can be obtained by determining all k-mers in the forward and in the backward strand for one
window. Only the lexicographically smallest k-mer in one window is saved, then the window is shifted by one and the
procedure is repeated until all windows have been processed. If two consecutive windows share the same minimiser, it is
stored only once.

The following figure shows three examples, where a) and b) differ only in their window size, while c) introduces
a gapped k-mer. All positions marked with a ’.’ are gaps and all k-mers in lower case come from the reverse strand.
As the example shows, k-mer size, shape and window size influence the total amount of minimisers.

\image html minimisers.png

### Non-lexicographical Minimisers

When sliding the window over the sequence, it might happen that consecutive minimisers differ only slightly.
For instance, when a minimiser starts with a repetition of A’s, it is highly likely that the minimiser of the next
window will start with a repetition of A’s too, which will then just be one A shorter. This may go on for multiple
window shifts, depending on how long the repetition is. Saving these only slightly different minimisers makes no sense
because they contain no new information about the underlying sequence.
Additionally, sequences with a repetition of A’s will be seen as more similar to each other than they actually are.
As [Marçais et al.](https://doi.org/10.1093/bioinformatics/btx235) have shown, randomizing the order of the k-mers
can solve this problem.

### Robust Winnowing

In case there are multiple minimal values within one window, the minimum and therefore the minimiser is ambiguous.
We choose the rightmost value as the minimiser of the window, and when shifting the window, the minimiser is only
changed if there appears a value that is strictly smaller than the current minimum. This approach is termed
*robust winnowing* by [Chirag et al.](https://doi.org/10.1093/bioinformatics/btaa435) and is
proven to work especially well on repeat regions.

# Usage in SeqAn3

The minimisers are implemented in `seqan3::views::minimiser_hash`. The function requires three arguments: The sequence,
the shape and the window size. The above-mentioned k-mer size is implicitly given by the size of the shape.
That is all the information you need to obtain minimisers with SeqAn3, so let's just dive right into the first
assignment!

\assignment{Assignment 1: Fun with minimisers I}
Task: Obtain the minimisers for "CCACGTCGACGGTT" with an ungapped shape of size 4 and a window size of 8.
\endassignment
\solution
\include doc/tutorial/06_minimisers/minimisers_solution1.cpp
\endsolution

If you have completed the assignment, you're probably wondering what those large numbers mean. As explained above,
lexicographical ordering is less than optimal. Therefore, SeqAn3 uses a seed to randomize the order. To do so, SeqAn3
simply XORs the hash values with a random seed (Default: 0x8F3F73B5CF1C9ADE). How would you get back the actual hash
values then?
Well, you just use XOR again!

\include doc/tutorial/06_minimisers/seed_example.cpp

From these hash values, you can obtain the sequence they represent by transforming the numbers to base 4. (For
example, 134 is "2012" in base four and therefore represents "GACG".)

Now take a closer look at the resulting minimiser sequences. Are they what you
expected? Probably not since they do not correspond to the minimsers computed in our original example. Can you figure
out why?

\hint
The example above is based on a lexicographical ordering, our first assignment is not.
\endhint

So, what would we need to change to achieve the results from our example? We need to use a different seed. SeqAn3
makes this simple by allowing us to add another input parameter of type `seqan3::seed`.

\assignment{Assignment 2: Fun with minimisers II}
Task: Reproduce the results from the example above by obtaining the minimisers for "CCACGTCGACGGTT".

To do so, you have to know which seed to use.

\hint
Because the randomization is based on an XOR, you have to use a value X which XOR-ed with a hash value results in the
same hash value.
\endhint

\hint
10001 XOR 0 = 10001
\endhint

\endassignment
\solution
\include doc/tutorial/06_minimisers/minimisers_solution2.cpp
\endsolution

### Ignoring the backward strand

So far, we have always considered the forward and the backward strand of a sequence, but there might be cases where
the backward strand should not be considered. If this is the desired behaviour then the `seqan3::views::minimiser` needs
to be used. Unlike the `seqan3::views::minimiser_hash`, `seqan3::views::minimiser` does not hash the values for you, so
you have to do this yourself. But fear not, `seqan3::views::kmer_hash` makes this really easy for you!

\snippet doc/tutorial/06_minimisers/minimisers_snippets.cpp minimiser


This syntax will result in minimisers with k-mer size 4 and a window-length of 8 (5 + 4 - 1). (So, to determine the
right window size you always have to determine the window size you would like to use in `seqan3::views::minimiser_hash`
and subtract it by the k-mer size and add 1, here: 8 - 4 + 1 = 5.)
Got it? Then let's make a simple test, what is the actual window size of the following:


```cpp
auto minimisers = text | seqan3::views::kmer_hash(seqan3::ungapped{4}) | seqan3::views::minimiser(1);
```
\hint
1 + 4 - 1 = 4. So, k-mer size and window size are equal.
\endhint

Wait a second, what does it mean, if k-mer size and window size are equal and we are not considering the backward
strand?

\hint
`seqan3::views::kmer_hash(seqan3::ungapped{4}) | seqan3::views::minimiser(1)` will return the same values as
`seqan3::views::kmer_hash`, because if k-mer size and window size are equal, there is only one comparison to determine
the minimal k-mer, the comparison between forward and backward strand. However, if the backward strand is not
considered, the minimisers are all k-mers from the given sequence.
\endhint

In order to ensure that this is the desired behaviour, using `seqan3::views::minimiser(1)` is prohibited. Please, use
`seqan3::views::kmer_hash` instead.

Last but not least, `seqan3::views::kmer_hash` and `seqan3::views::minimiser` do not have a seed parameter. So, in order
to obtain a random ordering, you have to XOR the view yourself. This can be done with the following command:

\snippet doc/tutorial/06_minimisers/minimisers_snippets.cpp minimiser_seed

\assignment{Assignment 3: Fun with minimisers III}
Task: Repeat assignment 2 but this time do not consider the backward strand.
\endassignment
\solution
\include doc/tutorial/06_minimisers/minimisers_solution3.cpp
\endsolution

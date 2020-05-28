# Minimiser {#tutorial_minimiser}

[TOC]

This tutorial introduces minimisers. Minimisers are a compact representation of a sequence that is closely related to,
but more efficient than a representation by k-mers. The minimisers are implemented as a view in SeqAn3. For more
information about views and how to implement your own view, please see  
\ref tutorial_ranges and \ref howto_write_a_view .

\tutorial_head{Easy, 20 min, , \ref tutorial_ranges\, \ref howto_write_a_view}

# Motivation

A common way to work with sequences is to obtain their k-mers, but even with a large size of k, the number of k-mers is
non-trivial, therefore storing k-mers is costly. At the same time, due to the overlap of consecutive k-mers, the
information they contain is highly redundant. That is why, minimiser come in handy. Minimisers are k-mers, which are
minimal in a window of a specified size.  Minimal could for example mean lexicographically smallest. By storing only
these minimal k-mers the storage cost is significantly reduced while maintaining a similar percentage of information.

# Minimiser Workflow

Because minimisers are minimal k-mers, they depend on a given k-mer size, a given shape (specifying which positions
should be considered) and a window size, which has to be greater or equal than the k-mer size. If all these values are
given, then the minimisers can be obtained by determining all k-mers in the forward and in the backward strand for one
window. Only the lexicographically smallest k-mer in one window is saved, then the window is shifted by one and the
procedure is repeated until all windows have been processed. If two consecutive windows share the same minimiser, it is
stored only once.

The following figure shows three examples, where a) and b) differ only in their window size, while c) introduces
a gapped k-mer. All positions marked with a ’.’ are gaps and all k-mers in lower case come from the reverse strand.
As the example shows, k-mer size, shape and window size influence the total amount of minimisers.

<img src="Minimiser.png"
     alt="Minimiser Example"
     style="width: 1000px; float: left; margin-right: 2px;" />


### Non-lexicographical Minimisers

When sliding the window over the sequence, it might happen that consecutive minimisers differ only slightly.
For instance, when a minimiser starts with a repetition of A’s, it is highly likely that the minimiser of the next
window will start with a repetition of A’s too, which will then just be one A shorter. This may go on for multiple
window shifts, depending on how long the repetition is. Saving these only slightly different minimiser makes no sense
because they contain no new information about the underlying sequence.
Additionally, sequences with a repetition of A’s will be seen as more similar to each other than they actually are.
As [Marçais et al.](https://doi.org/10.1093/bioinformatics/btx235) have shown, randomizing the order of the k-mers
can solve this problem.


# Usage in SeqAn3

The minimisers are implemented in `seqan3::views::minimiser_hash`. The function requires three arguments: The sequence,
the shape and the window size. The above mentioned k-mer size is implicitly given by the size of the shape.
That is all the information you need to obtain minimisers with SeqAn3, so let's just dive right into the first
assignment!

\assignment{Assignment 1: Fun with minimisers I}
Task: Obtain the minimisers for "CCACGTCGACGGTT" with an ungapped shape of size 4 and a window size of 8.
\endassignment
\solution
\include minimiser_solution1.cpp
\endsolution

If you have completed the assignment, you probably wonder what these large numbers mean. As explained above, the
lexicographical ordering is less optimal, therefore, SeqAn3 uses a seed to randomize the order. To do so, SeqAn3 simply
XORs the hash values with a random seed (Default: 0x8F3F73B5CF1C9ADE). How would you get back the actual hash values
then?
Well, you just use XOR again!

```cpp
uint64_t seed = 0x8F3F73B5CF1C9ADE;
auto hash_values = minimisers | std::views::transform([seed] (uint64_t i) {return i ^ seed;});
seqan3::debug_stream << hash_values << '\n'; // results in: [182, 216, 134]
```

From these hash values, you can obtain the sequence they are representing by transforming the numbers to base 4. (For
example, 182 is "2312" in base four and therefore represents "GTCG".)

Now take a closer look at the resulting minimiser sequences. Are they what you would
expect? Probably not, these are not the values presenting the minimisers shown in the example above. Why not?

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
\include minimiser_solution2.cpp
\endsolution

### Ignoring the backward strand

So far, we have always considered the forward and the backward strand of a sequence, but there might be cases where
the backward strand should not be considered. Therefore, there is an option to use `seqan3::views::minimiser_hash` with
the option `reverse==false`. What happens now, if we use `reverse==false` and the `shape.size()` equals the
`window_size`?

\hint
`seqan3::views::minimiser_hash` will return the same values as `seqan3::views::kmer_hash`, because if `shape.size()`
equals the `window_size` there is only one comparison to determine the minimal k-mer, the comparison between forward
and backward strand. Therefore, if the backward strand is not considered, the minimisers are all k-mers from the given
sequence.
\endhint

In order to ensure that this is the desired behaviour, using `seqan3::views::minimiser_hash` with `reverse==false`
and `shape.size()` equals `window_size` is prohibited. Please, use `seqan3::views::kmer_hash` instead.

\assignment{Assignment 3: Fun with minimisers III}
Task: Repeat assignment 2 but this time do not consider the backward strand.
\endassignment
\solution
\include minimiser_solution3.cpp
\endsolution

# First steps with SeqAn3 {#tutorial_first_example}

***Learning Objective:***

This tutorial walks you through small SeqAn3 programs. It is intended to give you a short overview
of what to expect in the other tutorials and how to use this documentation.

\tutorial_head{Easy, 30 min, \ref setup, [Ranges](https://github.com/seqan/seqan3/wiki/Ranges)\,
                                         [Concepts](https://en.cppreference.com/w/cpp/language/constraints)}

*Every page in the tutorials begins with this section. It is recommended that you do the "prerequisite tutorials"
before the current one. You should also have a look at the links provided in "recommended reading" and maybe keep
them open in separate tabs/windows as reference.*

***These tutorials try to briefly introduce C++ features not well known, however they do not teach programming in C++!
If you know how to program in another language, but are not familiar with C++ and/or the significant
changes in the language in recent years, we recommend the following resources:***

  * Bjarne Stroustrup: "A Tour of C++", Second Edition, 2018.

[TOC]

<br>

---

# Hello World!

Most good tutorials start with an easy *Hello World!* program. So have a look:

\snippet introduction_hello_world.cpp hello

Now you may ask, why we do not use std::cout or std::cerr for console output.
Actually, for the given text it does not make a difference since seqan3::debug_stream prints to std::cerr as well.
However, the debug stream provides convenient output for SeqAn's types as well as widely used data structures
(e.g. std::vector), which is especially helpful when you debug or develop your program
(that's where the name originates).

\assignment{Exercise: Debug stream}
Write a program that creates a std::vector of type `int` and initialise the vector with a few values.
Then print the vector with seqan3::debug_stream. Does your program also work with std::cerr?
\endassignment
\solution
\snippet introduction_debug_stream.cpp debug
\endsolution

\note
You may have spotted that the blue coloured keywords link you directly to the respective API documentation.
This is helpful if you need further information on a function, concept or class. We recommend you to open them
in separate browser tabs such that you can easily switch back to the tutorial.

<br>

---

# Read sequence files

SeqAn3 provides fast and easy access to biological file formats.
The following code example demonstrates the interface of seqan3::sequence_file_input.

\snippet introduction_align.cpp sequence_input_include
\snippet introduction_align.cpp sequence_input

Can you imagine anything easier? After you have initialised the instance with a filename,
you can simply query the data of interest using the get function. The available fields are
`SEQ`, `ID` and `QUAL` to retrieve sequences, ids and qualities, respectively. The latter is empty
unless you read FastQ files. The appropriate file format is detected by SeqAn from your file name suffix.

Here is the content of `seq.fasta`, so you can try it out!

~~~
>seq1
ACGTGATG
>seq2
AGTGATACT
~~~

<br>

---

# Align two sequences

We have two sequences from the file above now – so let us align them.
The pairwise sequence alignment is one of the core algorithms in SeqAn and used by several library components
and apps. It is strongly optimised for speed and parallel execution while providing exact results and a
generic interface.

\snippet introduction_align.cpp alignment_include
\snippet introduction_align.cpp alignment

The algorithm returns a range of result objects – which is the reason for the loop here (in this case the range
has length 1). Instead of passing a single pair of sequences, we could give a vector of sequence pairs to the
algorithm which then executes all alignments in parallel and stores the results in various seqan3::align_result
objects. The second argument to seqan3::align_pairwise is the *configuration* which allows you to specify
a lot of parameters for the alignment computation, for instance score functions, banded alignment and whether
you wish to compute a traceback or not. The configurations have their own namespace seqan3::align_cfg and can
be piped for building combinations. Check the alignment tutorial if you want to learn more.

\note
Every directory in the SeqAn sources contains an `all.hpp` file which includes all the functionality
of the respective (sub-) module. Thus, if you are working with several alignment methods, you can simply
include `seqan3/alignment/all.hpp` instead of naming each file individually.

<br>

---

# Indexed searching
At the end of this introduction tutorial, we want to present you another efficient algorithm in SeqAn:
The following code example shows you the steps for performing rapid searches in an indexed genome sequence.
As index we use the FM index here. After creation we want to query the positions of the motif *TAG* inside
a given sequence.

\snippet introduction_align.cpp index_search_include
\snippet introduction_align.cpp index_search

For more information on indexing and searching we refer you to the *Index Tutorial*. We hope you got a rough
overview on how the SeqAn3 interface has developed since version 2 and you are now free to continue with the
tutorials of your interest.

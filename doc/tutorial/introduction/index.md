# First steps with SeqAn {#tutorial_first_example}

***Learning Objective:***

This tutorial walks you through small SeqAn programs. It is intended to give you a short overview
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

# Hello World!

Most good tutorials start with an easy *Hello World!* program. So have a look:

\snippet introduction_hello_world.cpp hello

\note
This is a code snippet. You will see many code snippets in our documentation.
Most of them are compilable as-is, but some are only valid in their context,
e.g. they depend on other code snippets given before/after the current one or
other statements implied by the text. You can **copy'n'paste** freely from these examples,
this implies no copyright-obligations (however distributing SeqAn or an application
using it does, see [Copyright](https://docs.seqan.de/seqan/3-master-user/about_copyright.html) and [Citing](https://docs.seqan.de/seqan/3-master-user/about_citing.html)).

You may ask, why we do not use std::cout or std::cerr for console output.
Actually, for the given text it does not make a difference since seqan3::debug_stream prints to std::cerr as well.
However, the debug stream provides convenient output for SeqAn's types as well as widely used data structures
(e.g. std::vector), which is especially helpful when you debug or develop your program
(that's where the name originates).

\assignment{Assignment 1: Debug stream}
Write a program that creates a std::vector of type `int` and initialise the vector with a few values.
Then print the vector with seqan3::debug_stream. Does your program also work with std::cerr?
\endassignment
\solution
\snippet introduction_debug_stream.cpp debug
\endsolution

\note
This is an assignment with solution. You will find assignments in the tutorials to practise the discussed contents.
We believe that programming them will help you to memorise better and makes the tutorials more interesting and
interactive. The solutions provide the intended use; but often there are multiple ways to solve an assignment,
so don't worry too much if your solution is different from ours.

# Parse command line arguments

After we have seen the *Hello World!* program, we want to go a bit further and parse arguments from the command line.
The following snippet shows you how this is done in SeqAn. Here the program expects a string argument in the
program call and prints it to your terminal.

\snippet introduction_argument_parser.cpp argparse

Implementing a program with seqan3::argument_parser requires three steps:
1. Initialise the seqan3::argument_parser with your program's name and pass the `argc` and `argv` variables.
2. Register (positional) options in the parser object. In this way it knows which options to expect and
   it can generate the help page for your program. You will learn more about the option types in the *Argument Parser
   Tutorial*.
3. Run the parser. As it throws exceptions on wrong user behaviour, it should be surrounded with a try-catch block.

You will see that the entered text is now in the buffer variable `input`. The argument parser provides way more
functionality than we can show at this point, e.g. validation of arguments and different option types. We refer you
to the respective tutorial if you want to know more.

\note
You may have spotted that the blue coloured keywords link you directly to the respective **API documentation**.
This is helpful if you need further information on a function, concept or class. We recommend you to open them
in separate browser tabs such that you can easily switch back to the tutorial.

## Modules in SeqAn

You have just been introduced to one of the **Modules** of SeqAn, the *Argument Parser*.
Modules structure the SeqAn library into logical units, as there are for instance `alignment`, `alphabet`,
`argument_parser`, `io`, `search` and some more. See the *API Reference (Modules)* section in the
navigation column for a complete overview.

Some modules consist of submodules and the module structure is represented by the file hierarchy in the `include`
directory. Whenever you use functions of a module, make sure to `include` the correct header file.
Each directory in the SeqAn sources contains an `all.hpp` file which includes all the functionality
of the respective (sub-) module.
For small examples and quick prototyping, you can just include these `all.hpp`-headers.
However, for larger projects we recommend you include only the necessary headers, because this will reduce the
compile time measurably.

\note
If you remember the name of a function or class, but don't know which (sub-)module it belongs to,
you can enter it in the search bar (top-right).

# Read sequence files

Let's look at some functions of the IO module: SeqAn provides fast and easy access to biological file formats.
The following code example demonstrates the interface of seqan3::sequence_file_input.

\snippet introduction_file_input.cpp fileinput

Can you imagine anything easier? After you have initialised the instance with a filename,
you can simply step through the file in a for loop and retrieve the fields via
[structured bindings](https://en.cppreference.com/w/cpp/language/structured_binding).
The returned fields are `SEQ`, `ID` and `QUAL` to retrieve sequences, ids and qualities, respectively.
The latter is empty unless you read FastQ files. The appropriate file format is detected by SeqAn from
your filename's suffix.

Here is the content of `seq.fasta`, so you can try it out!

~~~
>seq1
ACGTGATG
>seq2
AGTGATACT
~~~

\assignment{Assignment 2: Read a FastA file}
Combine the code from above to read a FastA file and store its sequences in a std::vector of type seqan3::dna5_vector
(which is a common DNA sequence type in SeqAn). Use the argument parser for obtaining the filename as command line
argument to your program (e.g. call `./myprogram seq.fasta`).
\endassignment
\solution
\snippet introduction_read_fasta.cpp read
\endsolution

Note that the same code can also read FastQ files and the `qual` variable will not be empty then. If you like, try it!

\note
SeqAn3 uses `snake_case` for almost everything, also class names. Only C++ concepts are named using `CamelCase`.

# Align two sequences

We have two sequences from the file above now – so let us align them.
The pairwise sequence alignment is one of the core algorithms in SeqAn and used by several library components
and apps. It is strongly optimised for speed and parallel execution while providing exact results and a
generic interface.

\snippet introduction_align.cpp alignment

The algorithm returns a range of result objects – which is the reason for the loop here (in this case the range
has length 1). Instead of passing a single pair of sequences, we could give a vector of sequence pairs to the
algorithm which then executes all alignments in parallel and stores the results in various seqan3::alignment_result
objects. The second argument to seqan3::align_pairwise is the *configuration* which allows you to specify
a lot of parameters for the alignment computation, for instance score functions, banded alignment and whether
you wish to compute a traceback or not. The configurations have their own namespace seqan3::align_cfg and can
be combined via the logical OR operator (`|`) for building combinations. Check out the alignment tutorial if you want
to learn more.

\note
We encourage you to avoid declaring `using namespace seqan3;`. This has the additional benefit of easily distinguishing
between library features and standard C++. The only exception are string literals, where we often use
`using seqan3::operator""_dna4;` for convenience.

\note
We use a lot of Modern C++ in SeqAn so some things might look alien at first,
e.g. type templates are used like ordinary types in many situations (no `<>`).
We also always use `{}` to initialise objects and not `()` which is only used for function calls.
In general the style should be much easier for newcomers.

Now that you reached the end of this first tutorial, you know how SeqAn code looks like and you are able
to write some first code fragments. Let's go more into detail with the module-based tutorials!

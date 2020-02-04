# Parsing command line arguments with SeqAn {#tutorial_argument_parser}

<b>Learning Objective:</b> <br>
You will learn how to use the seqan3::argument_parser class to parse command line arguments. This tutorial is a walkthrough with links to the API documentation and is also meant as a source for copy-and-paste code.

\tutorial_head{Easy, 30-60 min, \ref setup, [POSIX conventions](https://www.math.uni-hamburg.de/doc/java/tutorial/essential/attributes/_posix.html)}

[TOC]

<br><br>

# Introduction

An easy and very flexible interface to a program is through the command line. This tutorial explains how to parse the command line using the SeqAn library’s seqan3::argument_parser class.

This class will give you the following functionality:

* Robust parsing of command line arguments.
* Simple validation of arguments (e.g. within a range or one of a list of allowed values).
* Automatically generated and nicely formatted help screens when your program is called with `--help`. You can also export this help to HTML and man pages.
* In the future, you are also able to automatically generate nodes for work flow engines such as KNIME or Galaxy.

## Command line argument terminology

Before we start, let's agree on some terminology. Consider the following command line call
```sh
mycomputer$ ./program1 -f -i 4 --long-id 6 file1.txt
```
The binary `program1` is called with several command line **arguments**. We call every single input an **argument** but differentiate their purpose into **options**, **positional_options**, **flags** or simply a **value** corresponding to one of the former. In our example above, `-f` is a **flag** which is never followed by a value, `-i` is an **option** with a short identifier (id) followed by its **value** `4`, `--long-id` is also an **option** with a long identifier followed by its **value** `6`, and **file1.txt** is a **positional_option**, because it is an option identified by its position instead of an identifier.

| Name                   | Purpose                                        | Example                 |
|:-----------------------|:-----------------------------------------------|:------------------------|
| **option**             | identify an argument by name (*id-value pair*) | `-i 5` or `--long-id 5` |
| **flag**               | boolean on/off flag (*id*)                     | `-f`                    |
| **positional option**  | identify an argument by position (*value*)     | `file1.txt`             |

Have a look at the [POSIX conventions](https://www.math.uni-hamburg.de/doc/java/tutorial/essential/attributes/_posix.html) for command line arguments if you want a detailed description on the requirements for the above. (Note: in the linked article the following holds: value="argument", option="option", flag = "option that does not require arguments", positional option ="non-option").

## A continuous example
We will get to know the wide functionality of the argument parser by writing a little application and extending it step by step.  Let's say we have a tab separated file `data.tsv` with information on the Game of Thrones Seasons ([by Wikipedia](https://en.wikipedia.org/wiki/List_of_Game_of_Thrones_episodes)):

\include doc/tutorial/argument_parser/data.tsv

We want to build an application that is able to read the file with or without a header line, select certain seasons and compute the average or median from the "Avg. U.S. viewers (millions)" of the selected seasons.

# The SeqAn argument parser class

Before we add any of the options, flags, and positional options, we will take a look at the seqan3::argument_parser class itself. It is constructed by giving a program's name and passing the parameters `argc` and `argv` from main. Note that no command line arguments have been parsed so far, but we can now add more information to the parser. After adding all desired information, the parsing is triggered by calling the seqan3::argument_parser::parse member function. Since the function throws in case any errors occur, we need to wrap it into a try-catch block. Here is a first working example:

\include doc/tutorial/argument_parser/basic_parser_setup.cpp

There are two types of exceptions: The seqan3::design_error which indicates that the parser setup was wrong (directed to the developer of the program, not the user!) and any other exception derived from seqan3::argument_parser_error, which detects corrupted user input. Additionally, there are special user requests that are handled by the argument parser by exiting the program via std::exit, e.g. calling `--help` that prints a help page screen.

## Design restrictions (seqan3::design_error)

The argument parser checks the following restrictions and throws a seqan3::design_error if they are not satisfied:

* **Long identifiers**: must be unique, more than one character long, may only contain alphanumeric characters, as well as `_`, `-`, or `@`, but never start with `-`.
* **Short identifiers**: must be unique and consist of only a single letter that is alphanumeric characters, `_` or `@`.
* either the short or long id may be empty but not both at the same time.
* Only the last positional option may be a list (see [lists](#section_list_positional_options)).
* The flag identifiers `-h`, `--help`, `--advanced-help`, `--advanced-help`, `--export-help`, `--version`, `--copyright` are predefined and cannot be specified manually or used otherwise.
* The seqan3::argument_parser::parse function may only be called once (per parser).

## Input restrictions

When calling the seqan3::argument_parser::parse function, the following potential user errors are caught (and handled by throwing a corresponding exception):

<table border="0">
<tr> <td> seqan3::unknown_option </td> <td> The option/flag identifier is not known to the parser. </td> </tr>
<tr> <td> seqan3::too_many_arguments </td> <td> More command line arguments than expected are given. </td> </tr>
<tr> <td> seqan3::too_few_arguments </td> <td> Less command line arguments than expected are given. </td> </tr>
<tr> <td> seqan3::required_option_missing </td> <td> A required option is not given (see [Required options](#section_required_option))</td> </tr>
<tr> <td> seqan3::user_input_error </td> <td> The given (positional) option value was invalid. </td> </tr>
<tr> <td> seqan3::validation_error </td> <td> (Positional-)Option validation failed (see [Validators](#section_validation)) </td> </tr>
</table>

## Special Requests (std::exit)

We denote "special requests" to command line input that does not aim to execute your program but rather display information about your program. Because on those request we never expect that the program is intended to run, we exit the program at the end of the seqan3::argument_parser::parse call via std::exit.

Currently we support the following *special requests*:

<table border="0">
<tr> <td> `-h/--help`    </td> <td> Prints the help page to the command line (`std::cout`) </td> </tr>
<tr> <td> `-hh/--advanced-help` </td> <td> Prints the advanced help page to the command line (`std::cout`) </td> </tr>
<tr> <td> `--export-help` </td> <td> Exports the help page in a different format (`std::cout`) </td> </tr>
<tr> <td> `--version`    </td> <td> Prints the version information to the command line (`std::cout`) </td> </tr>
<tr> <td> `--copyright`  </td> <td> Prints the copyright information to the command line (`std::cout`) </td> </tr>
</table>

\assignment{Assignment 1}
Copy the minimal working example into a cpp file in your working directory and compile it.
Play around with the binary, e.g. requesting special behaviour like printing the help page.
\endassignment

## Meta data {#section_meta_data}

Of course there is not much information to display yet, since we did not provide any. Let's improve this by modifying the seqan3::argument_parser::info member of our parser. The seqan3::argument_parser::info member is a struct of type seqan3::argument_parser_meta_data and contains the following members that can be customised:

- **app_name**, which is already set on construction (seqan3::argument_parser_meta_data::app_name)
- **author** (seqan3::argument_parser_meta_data::author)
- **citation** (seqan3::argument_parser_meta_data::citation)
- **date** (seqan3::argument_parser_meta_data::date)
- **description** (seqan3::argument_parser_meta_data::description)
- **email** (seqan3::argument_parser_meta_data::email)
- **examples** (seqan3::argument_parser_meta_data::examples)
- **long_copyright** (seqan3::argument_parser_meta_data::long_copyright)
- **man_page_title** (seqan3::argument_parser_meta_data::man_page_title)
- **short_copyright** (seqan3::argument_parser_meta_data::short_copyright)
- **short_description** (seqan3::argument_parser_meta_data::short_description)
- **synopsis** (seqan3::argument_parser_meta_data::synopsis)
- **url** (seqan3::argument_parser_meta_data::url)
- **version** (seqan3::argument_parser_meta_data::version)

\assignment{Assignment 2}

1. Extend the minimal example from assignment 1 by a function
   `void initialise_argument_parser(seqan3::argument_parser & parser)`.
2. Within this function, customise the parser with the following information:
   * Set the author to your favourite Game of Thrones character (Don't have one? Really? Take "Cersei").
   * Set the short description to "Aggregate average US. Game of Thrones viewers by season.".
   * Set the version to 1.0.0 .
   * Set some more, if you want to.
   Hint: Check out the API documentation for seqan3::argument_parser_meta_data and seqan3::argument_parser::info.
3. Try calling `--help` again and see the results.

\endassignment
\solution

\include doc/tutorial/argument_parser/solution1.cpp

\endsolution

# Adding options, flags and positional_options

Now that we're done with the meta information, we will learn how to add the actual functionality of options, flags and positional options. For each of these three there is a respective member function:

* seqan3::argument_parser::add_option
* seqan3::argument_parser::add_flag
* seqan3::argument_parser::add_positional_option

Each of the functions above take a variable by reference as the first parameter, which will directly store the corresponding parsed value from the command line. This has two advantages compared to other command line parsers: (1) There is no need for a getter function after parsing and (2) the type is automatically deduced (e.g. with boost::program_options you would need to access `parser["file_path"].as<std::filesystem::path>()` afterwards).

The seqan3::argument_parser::add_flag only allows boolean variables while seqan3::argument_parser::add_option and seqan3::argument_parser::add_positional_option allow **any type that is convertible from a std::string via std::from_chars** or a container of the former (see \ref tutorial_argument_parser_list_options). Besides accepting generic types, the parser will **automatically check if the given command line argument can be converted into the desired type** and otherwise throw a seqan3::type_conversion_error exception.

So how does this look like? The following code snippet adds a positional option to `parser`.

\snippet doc/tutorial/argument_parser/small_snippets.cpp add_positional_option

Additionally to the variable that will store the value, you need to pass a description. This description will help users of your application to understand how the option is affecting your program.

\note As the name suggest, positional options are identified by their position. In SeqAn, the first `add_positional_option()` will be linked to the first command line argument that is neither an option-value pair nor a flag. So the order of initialising your parser determines the order of assigning command line arguments to the respective variables.
We personally recommend to always use regular options (id-value pairs) because they are more expressive and it is easier to spot errors.

You can add an option like this:

\snippet doc/tutorial/argument_parser/small_snippets.cpp add_option

Additionally to the variable that will store the value and the description, you need to specify a short and long identifier. The example above will recognize an option `-n` or `--my-number` given on the command line and expect it to be followed by a value separated only by `=` or space or by nothing at all.

Finally, you can add a flag with the following call:

\snippet doc/tutorial/argument_parser/small_snippets.cpp add_flag

Note that you can omit either the short identifier by passing <code>'\0'</code> or the long identifier by passing `""` but you can never omit both at the same time.

## Default values

With the current design, every option/flag/positional automatically has a **default value** which simply is the value with which you initialise the corresponding variable that is passed as the first parameter. Yes it is that easy, just make sure to always initialise your variables properly.

\assignment{Assignment 3}
Getting back to our example application, let's extend our code from Assignment 2 by the following:

As a *best practice* recommendation for handling multiple options/flags/positionals, you should store the variables in a struct and pass this struct to your parser initialisation function. You can use the following code that introduces a `cmd_arguments` struct storing all relevant command line arguments. Furthermore, it provides you with a small function `run_program()` that reads in the data file and aggregates the data by the given information. You don't need to look at the code of `run_program()`, it is only so that we have a working program.

\htmlonly <details> <summary>Copy and paste this code into the beginning of your application: </summary>  \endhtmlonly
\snippet doc/tutorial/argument_parser/solution3.cpp program
\htmlonly </details> \endhtmlonly

Your task is now to extend the initialisation function by the following:

1. Extend your `initialise_argument_parser` function by a parameter that takes a `cmd_arguments` object and adapt the function call in your main function to pass on `args`;
2. Set the default value of `aggregate_by` to `"mean"`.

You can now use the variables from `args` to add the following inside of the `initialise_argument_parser` function:

3. Add a positional option to the parser that sets the variable `file_path` so our program knows the location of the data file to read in.
4. Add an option `-y/--year` that sets the variable `year`, which will enable our program to filter the data by only including a season if it got released after the value `year`.
5. Add an option `-a/--aggregate-by` that sets the variable `aggregate_by`, which will enable our program choose between aggregating by mean or median.
6. Add a flag `-H/--header-is-set` that sets the variable `header_is_set`, which lets the program know whether it should ignore the first line in the file.

Take a look at the help page again after you've done all of the above. You will notice that your options have been automatically included. **Copy and paste the example data file from the introduction** and check if your options are set correctly by trying the following few calls:

```console
./game_of_parsing -H -y 2014 data.tsv
7.9175
```

```console
./game_of_parsing -H -y 2010 --aggregate-by median data.tsv
6.84
```
\endassignment
\solution
\snippet doc/tutorial/argument_parser/solution3.cpp solution

\htmlonly <div class=\"assignment\"> <details><summary>In case you are stuck, the complete program now looks like this:</summary> \endhtmlonly
\include doc/tutorial/argument_parser/solution3.cpp
\htmlonly </details> </div> \endhtmlonly

\endsolution

# List options {#tutorial_argument_parser_list_options}

In some use cases you may want to allow the user to specify an option multiple times and store the values in a list. With the seqan3::argument_parser this behaviour can be achieved simply by choosing your input variable to be of a container type (e.g. std::vector). The parser registers the container type through seqan3::container and will adapt the parsing of command line arguments accordingly.

Example:

\snippet doc/tutorial/argument_parser/small_snippets.cpp option_list

Adding this option to a parser will allow you to call the program like this:

```console
./some_program -n Jon -n Arya -n Ned
```

The vector `list_variable` will then contain all three names `["Jon", "Arya", "Ned"]`.

## List positional options? {#section_list_positional_options}

An arbitrary positional option cannot be a list because of the ambiguity of which value belongs to which positional option. We do allow the very last option to be a list for convenience though. Note that if you try to add a positional list option which is not the last positional option, a seqan3::design_error will be thrown.

Example:

\snippet doc/tutorial/argument_parser/small_snippets.cpp positional_option_list

Adding these positional options to a parser will allow you to call the program like this:

```console
./some_program Stark Jon Arya Ned
```

The first `variable` will be filled with the value `Stark` while the vector `list_variable` will then contain the three names `["Jon", "Arya", "Ned"]`.

\assignment{Assignment 4}

We extend the solution from assignment 3:

1. Remove the option `-y/--year`, since we want to keep it simple and only aggregate by season now.

2. Add a variable `seasons` of type `std::vector<uint8_t>` to the struct `cmd_arguments`.

3. Add a list option `-s/--season` that will fill the variable `seasons` which lets the user specify which seasons to aggregate instead of the year.

4. [BONUS] If you have some spare time, try to adjust the program code to aggregate by season. Hint: Use std::find.
\htmlonly <details> <summary>Otherwise just replace the while loop with the following: </summary>  \endhtmlonly
\snippet doc/tutorial/argument_parser/solution4.cpp altered_while
\htmlonly </details> \endhtmlonly

Take a look at the help page again after you've done all of the above. You will notice that your option `-s/--season` even tells you that it is of type `List of unsigned 8 bit integer's`. Check if your options are set correctly by trying the following few calls:

```console
./game_of_parsing -H -s 2 -s 4 data.tsv
5.32
```

```console
./game_of_parsing -H -s 1 --season 3 -s 7 --aggregate-by median data.tsv
4.97
```
\endassignment
\solution
\snippet doc/tutorial/argument_parser/solution4.cpp solution

\htmlonly <div class=\"assignment\"> <details><summary>In case you are stuck, the complete program now looks like this:</summary> \endhtmlonly
\include doc/tutorial/argument_parser/solution4.cpp
\htmlonly </details> </div> \endhtmlonly

\endsolution

# Setting options as required, advanced or hidden

## Required options {#section_required_option}

There is a flaw in the example application we have programmed in assignment 4, did you notice? You can make it misbehave by not giving it any option `-s` (which is technically correct for the seqan3::argument_parser because a list may be empty). You could of course handle this in the program itself by checking whether the vector `seasons` is empty, but since supplying no season is not expected we can force the user to supply the option at least once by **declaring an option as required**.

For this purpose we need to use the seqan3::option_spec enum interface that is accepted as an additional argument by all of the `add_[positional_option/option/flag]` calls:

\snippet doc/tutorial/argument_parser/small_snippets.cpp required_option

If the user does **not** supply the required option via the command line, he will now get the following error:

```console
./example_program --some-other-option
Option -n/--name is required but not set.
```

\note Positional options are always required!

## Advanced and hidden options

Additionally to the *required* tag, there is also the possibility of **declaring an option as advanced or hidden**.

Set an option/flag to advanced, if you do not want the option to be displayed in the normal help page (`-h/--help`). Instead, the advanced options are only displayed when calling `-hh/--advanced-help`. This can be helpful, if you want to avoid to bloat your help page with too much information for inexperienced users of your application, but still provide thorough information on demand.

Set an option/flag to hidden, if you want to completely hide it from the user. It will neither appear on the help page nor in any export format. For example, this might be useful for debugging reasons.

Summary:

| Tag        | Description                                                    |
|:-----------|:---------------------------------------------------------------|
| DEFAULT    | The default tag with no special behaviour.                     |
| REQUIRED   | Required options will cause an error if not provided.          |
| ADVANCED   | Advanced options are only displayed wit `-hh/--advanced-help`. |
| HIDDEN     | Hidden options are never displayed when exported.              |

\assignment{Assignment 5}

Extend the solution from assignment 4 by declaring the `-s/--season` option as required.

Check if your options are set correctly by trying the following call:

```console
./game_of_parsing -H --aggregate-by median data.tsv
[Winter has come] Option -s/--season is required but not set.
```
\endassignment
\solution
\snippet doc/tutorial/argument_parser/solution5.cpp solution

\htmlonly <div class=\"assignment\"> <details><summary>In case you are stuck, the complete program now looks like this:</summary> \endhtmlonly
\include doc/tutorial/argument_parser/solution5.cpp
\htmlonly </details> </div> \endhtmlonly
\endsolution

# Validation of (positional) option values {#section_validation}

Our applications often do not allow just any value to be passed as input arguments and if we do not check for them, the program may run into undefined behaviour. The best way to carefully restrict user input is to directly check the input when parsing the command line. The seqan3::argument_parser provides **validators** for a given (positional) option.

A *validator* is a [functor](https://stackoverflow.com/questions/356950/what-are-c-functors-and-their-uses) that is called within the argument parser after retrieving and converting a command line argument. We provide several validators, which we hope cover most of the use cases, but you can always create your own validator (see section [Create your own validator](#section_create_your_own_validator)).

\attention You can pass a validator to the seqan3::argument_parser::add_option function only after passing the seqan3::option_spec parameter. Pass the seqan3::option_spec::DEFAULT tag, if there are no further restrictions on your option.

## SeqAn validators

The following validators are provided in the SeqAn library and can be included with the following header:

\snippet doc/tutorial/argument_parser/small_snippets.cpp validator_include

All the validators below work on single values or a container of values. In case the variable is a container, the validator is called **on each element** separately.

\note If the validators below do not suit your needs, you can always create your own validator. See the concept tutorial for an example of how to create your own validator.

### The seqan3::arithmetic_range_validator

On construction, this validator receives a maximum and a minimum number.
The validator throws a seqan3::validation_error exception whenever a given value does not lie inside the given min/max range.

\snippet test/snippet/argument_parser/validators_1.cpp validator_call

Our application has a another flaw that you might have noticed by now: If you supply a season that is not in the data file, the program will again misbehave. Instead of fixing the program, let's restrict the user input accordingly.

\assignment{Assignment 6}
Add a seqan3::arithmetic_range_validator to the `-s/--season` option that sets the range to `[1,7]`.
\endassignment
\solution
\snippet doc/tutorial/argument_parser/solution6.cpp arithmetic_range_validator
\endsolution

### The seqan3::value_list_validator

On construction, the validator receives a list (vector) of valid values.
The validator throws a seqan3::validation_error exception whenever a given value is not in the given list.

\snippet test/snippet/argument_parser/validators_2.cpp validator_call

\assignment{Assignment 7}
Add a seqan3::value_list_validator to the `-a/--aggregate-by` option that sets the list of valid values to `["median", "mean"]`.
\endassignment
\solution
\snippet doc/tutorial/argument_parser/solution6.cpp value_list_validator
\endsolution

### The file validator

SeqAn offers two file validator types: the seqan3::input_file_validator and the seqan3::output_file_validator.
On construction, the validator receives a list (vector) of valid file extensions that are tested against the extension
of the parsed option value.
The validator throws a seqan3::validation_error exception whenever a given filename's extension is not in the
given list of valid extensions. In addition, the seqan3::input_file_validator checks if the file exists, is a regular
file and is readable.
The seqan3::output_file_validator on the other hand ensures that the output does not already exist (in order to prevent
overwriting an already existing file) and that it can be created.

\note If you want to allow any extension just use a default constructed file validator.

Using the seqan3::input_file_validator:

\snippet test/snippet/argument_parser/validators_input_file.cpp validator_call

Using the seqan3::output_file_validator:

\snippet test/snippet/argument_parser/validators_output_file.cpp validator_call

### The directory validator

In addition to the file validator types, SeqAn offers directory validator types. These are useful if one needs
to provide an input directory (using the seqan3::input_directory_validator) or output directory
(using the seqan3::output_directory_validator) where multiple files need to be read from or written to.
The seqan3::input_directory_validator checks whether the specified path is a directory and is readable.
Similarly, the seqan3::output_directory_validator checks whether the specified directory is writable and can be created,
if it does not already exists.
If the tests fail, a seqan3::validation_error exception will be thrown. Also, if something unexpected with the
filesystem happens, a std::filesystem_error will be thrown.

Using the seqan3::input_directory_validator:

\snippet test/snippet/argument_parser/validators_input_directory.cpp validator_call

Using the seqan3::output_directory_validator:

\snippet test/snippet/argument_parser/validators_output_directory.cpp validator_call

\assignment{Assignment 8}
Add a validator to the first positional option that expects a file formatted with tab separated values.
Store the result in `file_path`.
\endassignment
\solution
\snippet doc/tutorial/argument_parser/small_snippets.cpp input_file_validator
\endsolution

### The seqan3::regex_validator

On construction, the validator receives a pattern for a regular expression.
The pattern variable will be used for constructing an std::regex and the validator will call std::regex_match on the command line argument.

Note that a regex_match will only return true if the string matches the pattern completely (in contrast to regex_search which also matches substrings). The validator throws a seqan3::validation_error exception whenever a given parameter does not match the given regular expression.

\snippet test/snippet/argument_parser/validators_4.cpp validator_call

## Chaining validators

You can also chain validators using the pipe operator (`|`). The pipe operator is the AND operation for two validators, which means that a value must pass both validators in order to be accepted by the combined validator.

For example, you may want a file name that only accepts absolute paths, but also must have one out of a list of given file extensions.
For this purpose you can chain a seqan3::regex_validator to a seqan3::input_file_validator:

\snippet test/snippet/argument_parser/validators_chaining.cpp validator_call

You can chain as many validators as you want, they will be evaluated one after the other from left to right (first to last).

\assignment{Assignment 9}
Add a seqan3::regex_validator to the first positional option that expects the `file_path` by chaining it to the already present seqan3::input_file_validator. The parsed file name should have a suffix called `seasons`.
\endassignment
\solution
\snippet doc/tutorial/argument_parser/solution6.cpp file_validator
\endsolution

# Full solution

The following solution shows the complete code including all the little assignments of this tutorial
that can serve as a copy'n'paste source for your application.

\solution
\include doc/tutorial/argument_parser/solution6.cpp
\endsolution

# Subcommand argument parsing

Many applications provide several sub programs, e.g. `git` comes with many functionalities like `git push`,
`git pull`, `git checkout`, etc. each having their own help page.
If you are interested in how this subcommand parsing can be done with the seqan3::argument_parser,
take a look at our \link subcommand_arg_parse HowTo\endlink.

# Update Notifications

When you run a SeqAn-based application for the first time, you will likely be asked about "update notifications".
This is a feature that helps inform users about updates
and helps the SeqAn project get a rough estimate on which SeqAn-based apps are popular.

See the API documentation of seqan3::argument_parser for information on how to configure (or turn off) this feature.
See our [wiki entry](https://github.com/seqan/seqan3/wiki/Update-Notifications) for more information on how it works and our privacy policy.

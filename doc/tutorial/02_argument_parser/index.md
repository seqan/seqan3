# Parsing command line arguments with Sharg {#tutorial_argument_parser}

<!-- SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
     SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
     SPDX-License-Identifier: CC-BY-4.0
-->

We have separated the feature of parsing command line arguments to its own project:

### [![sharg_logo][sharg_logo_link]][sharg_link] The Sharg Parser

* Github Repository: https://github.com/seqan/sharg-parser
* API Documentation: https://docs.seqan.de/sharg.html
* Tutorials: https://docs.seqan.de/sharg/main_user/usergroup1.html

<!-- Use the Sharg logo from permalink. -->
[sharg_logo_link]: https://raw.githubusercontent.com/seqan/sharg-parser/1.0.0/test/documentation/sharg_logo.svg "Open Github"
<!-- Link the logo to the documentation website. -->
[sharg_link]: https://github.com/seqan/sharg-parser

---

## Sharg & SeqAn

You can easily setup Sharg parallel to SeqAn as we use the exact same infrastructure.

If you have completed the \ref setup, do the following to also include the Sharg parser:

1. In the `tutorial` directory, clone the Sharg parser
  ```
  git clone https://github.com/seqan/sharg-parser.git
  ```
  Your directory structure now looks like this:
  ```
  tutorial
  ├── build
  ├── seqan3
  └── sharg-parser
  └── source
  ```

2. Adapt your CMake script:

<!-- Parsing the snippet like this to avoid verbatim includes of the snippet identifiers if we used nested snippets. -->
<!-- Snippet start -->
\dontinclude test/external_project/seqan3_setup_tutorial_with_sharg/CMakeLists.txt
\skipline cmake_minimum_required
\until target_link_libraries
<!-- Snippet end -->

Done!

**Now you can do the [basic tutorial of the Sharg parser](https://docs.seqan.de/sharg/main_user/tutorial_parser.html)
to learn how to conveniently access command line arguments.**

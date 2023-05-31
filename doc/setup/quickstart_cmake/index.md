# Quick Setup (using CMake) {#setup}

<b>Learning Objective:</b><br>
In this short guide you will learn how to set up SeqAn and how to compile a small example to test whether everything
works.

\tutorial_head{Easy, 30 Minutes, No prerequisites, }

[TOC]

<br>

# Software
Requirements:
  - gcc >= 11
  - cmake >= 3.4
  - git

## Installing GCC

SeqAn requires GCC >= 11. Current versions of LLVM/Clang and VisualStudio/MSVC are **not yet supported**.
We will briefly explain how to install GCC-11 (or the latest GCC if such an option is available) on some popular
operating systems. We recommend using the latest version of GCC available. For more information, refer to your
operating system's documentation.

### Linux
<div class="tabbed">
- <b class="tab-title">Ubuntu >= 22.04</b>
```bash
# Installs default compiler version (gcc-11 for Ubuntu 22.04).
sudo apt install g++
# To install gcc-12, follow instructions for Ubuntu < 22.04.
```
- <b class="tab-title">Ubuntu < 22.04</b>
```bash
sudo add-apt-repository --no-update --yes ppa:ubuntu-toolchain-r/ppa
sudo add-apt-repository --no-update --yes ppa:ubuntu-toolchain-r/test
sudo apt-get update
sudo apt install g++-11
```
- <b class="tab-title">Using conda</b>
To avoid interference with system packages, we recommend creating a new environment when using conda.
```bash
conda create -n conda_gcc_env -c conda-forge gcc_linux-64
conda activate conda_gcc_env
```
This will put GCC in a separate environment `conda_gcc_env` which can be activated via `conda activate conda_gcc_env`
and deactivated via `conda deactivate`.

</div>

### macOS
<div class="tabbed">
- <b class="tab-title">Using Homebrew</b>
```bash
brew install gcc@11
```
- <b class="tab-title">Using Macports</b>
```bash
sudo port install gcc11
```

</div>

### Windows
The Windows Subsystem for Linux offers an easy way to run a Linux distribution under Windows.
Follow [Microsoft's setup guide](https://docs.microsoft.com/en-us/windows/wsl/about) to install WSL and then follow
the steps listed for Linux-based systems.

### Browser
<div class="tabbed">
- <b class="tab-title">Using gitpod.io</b>
[gitpod.io](https://gitpod.io) allows you to edit, compile and run code from within your browser. The free version includes 50
hours of use per month, which is plenty for our tutorials. A GitHub account is required.
[Click here](https://gitpod.io/#https://github.com/seqan/seqan3/) to open SeqAn3 in gitpod.
- <b class="tab-title">Using GitHub Codespaces</b>
[GitHub Codespaces](https://github.com/codespaces) offer a service similar to gitpod, including a free monthly quota.
[Click here](https://codespaces.new/seqan/seqan3) to open SeqAn3 in Codespaces.
Please note that you may have to manually clone the submodules by running `git submodule update --init`.

</div>
<br>
\attention After installing, `g++ --version` should print the desired GCC version.
           If not, you may have to use, for example, `g++-11 --version` or even specify the full path to your compiler.

Similarly, you may need to install CMake and git, e.g. `apt install cmake git`.

# Directory Structure
In this section we will use the `tree` command to show the directory structure. This program may not be installed
on your system. If so, you may wish to install it or verify the directory structure in other ways, e.g. by using
`ls -l`.

For this project, we recommend following directory layout:

```
tutorial
├── source
├── build
└── seqan3
```

To set these directories up you can follow this script (note the <b>\--recurse-submodules</b> when cloning SeqAn3):
```bash
mkdir tutorial
cd tutorial
mkdir build
mkdir source
git clone --recurse-submodules https://github.com/seqan/seqan3.git
```

The output of the command `tree -L 2` should now look like this:
```
.
├── build
├── seqan3
│   ├── CHANGELOG.md
│   ├── CMakeLists.txt
│   ├── ...
│   └── test
└── source

8 directories, 6 files
```

# Compiling and Running

To test whether everything works, we will now compile and run a small example.

First we create the file `hello_world.cpp` in the `source` directory with the following contents:

\include test/external_project/src/hello_world.cpp

To compile it, we first create a `CMakeLists.txt` file in the `source` directory:
<!-- Parsing the snippet like this to avoid verbatim includes of the snippet identifiers if we used nested snippets. -->
<!-- Snippet start -->
\dontinclude test/external_project/seqan3_setup_tutorial/CMakeLists.txt
\skipline cmake_minimum_required
\until target_link_libraries
<!-- Snippet end -->

The directories should now look like this:

```
.
├── build
├── seqan3
│   ├── CHANGELOG.md
│   ├── CMakeLists.txt
│   ├── ...
│   └── test
└── source
    ├── CMakeLists.txt
    └── hello_world.cpp
```

Now we can switch to the directory `build` and run:

```bash
cmake -DCMAKE_BUILD_TYPE=Release ../source
make
./hello_world
```

The output should be `Hello World!`. Note that the build type is specified with `-DCMAKE_BUILD_TYPE=Release`.
Specifying `Release` enables an optimized build where no debug information is available. Release mode is therefore
suitable for the end user. Programs built using `-DCMAKE_BUILD_TYPE=Debug` will run slower, but also make the detection
of errors easier. `Debug` is suitable for contributors, and we recommend using it while working with our
[Tutorials](usergroup1.html).

\anchor remark_cmake_cxx_compiler
\remark Depending on the standard C++ compiler on your system, you may need to specify the compiler via
`-DCMAKE_CXX_COMPILER=`, for example:
```bash
cmake -DCMAKE_CXX_COMPILER=/path/to/executable/g++-11 ../source
```

# Adding a new source file to your project

If you create a new `cpp` file and want to compile it, you need to add another `add_executable` and
`target_link_libraries` directive to you `CMakeLists.txt`.
For example, after adding `another_program.cpp` your `CMakeLists.txt` may look like this:
\snippet test/external_project/seqan3_setup_tutorial/CMakeLists.txt adding_files

# Including SeqAn3 as external project

```cmake
cmake_minimum_required (VERSION 3.14)

project (my_app LANGUAGES CXX VERSION 1.0.0)

set (seqan3_git_tag "#.#.#") # adapt as needed, e.g. "3.2.0" or "master"

message (STATUS "Fetching SeqAn3 ${seqan3_git_tag}:")

include (FetchContent)
FetchContent_Declare (
    seqan3_fetch_content
    GIT_REPOSITORY "https://github.com/seqan/seqan3.git"
    GIT_TAG "${seqan3_git_tag}"
)

# Download and make SeqAn3 available.
FetchContent_MakeAvailable (seqan3_fetch_content)

add_executable (my_app my_app.cpp)

# Set up everything needed to use SeqAn3 with my_app:
target_link_libraries (my_app PUBLIC seqan3::seqan3)
```

# Encountered issues

* **Using conda's gcc package:** ``/usr/lib/x86_64-linux-gnu/libstdc++.so.6: version 'CXXABI_1.3.11' not found``<br>
  Try setting `LD_LIBRARY_PATH`:
  ```bash
  export LD_LIBRARY_PATH=<conda_install_path>/envs/conda_gcc_env/lib/
  ```
  where `<conda_install_path>` must be replaced by the path yo your conda installation.<br>
  Usually this corresponds to the path printed by `conda info --base` and may look similar to `/home/user/miniconda3/`.

* **Assembler not found:** `... could not understand flag m ...`<br>
  Try adding `/usr/bin` to your `PATH`:
  ```bash
  export PATH=/usr/bin:$PATH
  ```
  and run `cmake` again.

* **SDSL library not found**: `The SDSL library is required, but wasn't found.`<br>
  The repository was not cloned correctly. This can be verified by checking whether
  `submodules/sdsl-lite/include/version.hpp` exists.
  If it does not, try running `git submodule update --init` within the seqan3 directory.

* **Incorrect compiler**: `Your compiler is not supported.` or `Only GCC is supported.`<br>
  The incorrect compiler is used (e.g., Apple Clang instead of GCC). Be sure to set `-DCMAKE_CXX_COMPILER=`. For an
  example, see \ref remark_cmake_cxx_compiler "this remark".

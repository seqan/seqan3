# Quick Setup (using CMake) {#setup}

<b>Learning Objective:</b><br>
In this short guide you will learn how to set up SeqAn and how to compile a small example to test whether everything
works.

\tutorial_head{Easy, 30 Minutes, No prerequisites, }

[TOC]

<br>

# Software
Requirements:
  - gcc >= 7
  - cmake >= 3.4
  - git

SeqAn3 requires a compiler with full C++20 support *or* GCC >= 7. Current versions of LLVM/Clang and VisualStudio/MSVC are **not supported**.
We will briefly explain how to install GCC-7 on some popular operating systems, but we recommend using the latest version of GCC available. For more information refer to your operating system's documentation.

**Ubuntu >= 18.04**
```
sudo apt install g++
```
**Ubuntu < 18.04**
```
sudo add-apt-repository ppa:ubuntu-toolchain-r/test
sudo apt update
sudo apt install g++-7
```
**MacOS** using [Homebrew](https://brew.sh/)
```
brew install gcc@7
```

**MacOS** using [MacPorts](https://www.macports.org/)
```
sudo port install gcc-7
```

**Linux** using [conda](https://conda.io)
```
conda create -n gcc7 -c quantstack gcc-7 libgcc-7
conda activate gcc7
```
This will put GCC-7 in a separate environment called `gcc7` which can be activated via `conda activate gcc7` and deactivated via `conda deactivate gcc7`.

\note \htmlonly <div class=\"assignment\"> <details><summary><b>Known Issue:</b></summary> \endhtmlonly If you encounter the error <code>/usr/lib/x86_64-linux-gnu/libstdc++.so.6: version 'CXXABI_1.3.11' not found</code>, you have to set the LD_LIBRARY_PATH:
```
export LD_LIBRARY_PATH=/home/user/miniconda3/envs/gcc7/lib/
```
where `/home/user/miniconda3/` is the path to your conda installation. \htmlonly </details> </div> \endhtmlonly

\attention After installing, `g++ --version` should print the desired version. If not, you may have to use, for example, `g++-7 --version` or even specify the full path to your executable.

Similarly you need to install cmake and git, e.g. `apt install cmake git`.

# Directory Structure
For this project, we recommend following directory layout:

```
tutorial
├── source
├── build
└── seqan3
```

To set these directories up you can follow this script (note the <b>\--recurse-submodules</b> when cloning SeqAn3):
```
mkdir tutorial
cd tutorial
mkdir build
mkdir source
git clone --recurse-submodules https://github.com/seqan/seqan3.git
```

The output of the command ``` tree -L 2 ``` should now look like this:
```
.
├── build
├── seqan3
│   ├── build_system
│   ├── CHANGELOG.md
│   ├── CMakeLists.txt
│   ├── CODE_OF_CONDUCT.md
│   ├── CONTRIBUTING.md
│   ├── doc
│   ├── include
│   ├── LICENSE.md
│   ├── README.md
│   ├── submodules
│   └── test
└── source

8 directories, 6 files
```

# Compiling and Running

To test whether everything works, we will now compile and run a small example.

First we create the file `hello_world.cpp` in the `source` directory with the following contents:

\include hello_world.cpp

To compile it we first create a `CMakeLists.txt` file in the `source` directory:

```cmake
cmake_minimum_required (VERSION 3.4)
project (seqan3_tutorial CXX)

find_package (SeqAn3 3.0.0 REQUIRED HINTS "${CMAKE_SOURCE_DIR}/../seqan3/build_system")

add_executable (hello_world hello_world.cpp)

target_link_libraries (hello_world seqan3::seqan3)
```

The directories should now look like this:

```
tutorial
├── source
    ├── CMakeLists.txt
    └── hello_world.cpp
├── build
└── seqan3
    ├── CMakeLists.txt
    ├── LICENSE
    ...
```

Now we can switch to the directory `build` and run:

```bash
cmake ../source
make
./hello_world
```

The output should be `Hello world`.

\remark Depending on the standard C++ on your system, you may need to specify the compiler via `-DCMAKE_CXX_COMPILER=`, for example:
```bash
cmake ../source -DCMAKE_CXX_COMPILER=/path/to/executable/g++-7
```

\note In some situations it can happen that the correct assembler is not found.
You will see an error during the cmake configuration that says something like: `... could not understand flag m ...`.
In this case you can try to export the Path:
```
export PATH=/util/bin:$PATH
```
and try running cmake again.

# Adding a new source file to your project

If you create a new `cpp` file and want to compile it, you need to add another `add_executable` and
`target_link_libraries` directive to you `CMakeLists.txt`.
For example, after adding `another_program.cpp` your `CMakeLists.txt` may look like this:
```cmake
cmake_minimum_required (VERSION 3.4)
project (seqan3_tutorial CXX)

find_package (SeqAn3 3.0.0 REQUIRED HINTS "${CMAKE_SOURCE_DIR}/../seqan3/build_system")

add_executable (hello_world hello_world.cpp)
add_executable (another_program another_program.cpp)

target_link_libraries (hello_world seqan3::seqan3)
target_link_libraries (another_program seqan3::seqan3)
```

![Zenodo Badge](https://zenodo.org/badge/DOI/10.5281/zenodo.3755142.svg)

Fill-Spill-Merge
===========================================

**Title of Manuscript**:
 * Computing water flow through complex landscapes, Part 3: Fill-Spill-Merge: Flow routing in depression hierarchies (doi: [10.5194/esurf-9-105-2021](https://dx.doi.org/10.5194/esurf-9-105-2021))

**Previous Manuscripts**:
 * Computing water flow through complex landscapes, Part 2: Finding hierarchies in depressions and morphological segmentations (doi: [10.5194/esurf-2019-34](https://doi.org/10.5194/esurf-2019-34))
 * Computing water flow through complex landscapes – Part 1: Incorporating depressions in flow routing using FlowFill (doi: [10.5194/esurf-7-737-2019](https://doi.org/10.5194/esurf-7-737-2019))

**Authors**: Richard Barnes, Kerry Callaghan, Andrew Wickert

**Corresponding Author**: Richard Barnes (richard.barnes@berkeley.edu)

**Code Repositories**
 * [Author's GitHub Repository](https://github.com/r-barnes/Barnes2020-FillSpillMerge)

This repository contains a reference implementation of the algorithms presented
in the manuscript above, along with information on acquiring the various
datasets used, and code to perform correctness tests.



Abstract
--------

Depressions---inwardly-draining regions---are common to many landscapes,
including those shaped by glaciation, compressional and/or extensional
tectonics, and cratering. When there is sufficient moisture, depressions take
the form of lakes and wetlands; otherwise, they may be dry. Hydrological flow
models used in geomorphology, hydrology, planetary science, soil and water
conservation, and other fields often eliminate depressions through filling or
breaching; however, this can produce unrealistic results. Models that retain
depressions, on the other hand, are often undesirably expensive to run. In
previous work we began to address this by developing a depression hierarchy data
structure to capture the full topographic complexity of depressions in a region.
Here, we extend this work by presenting a Fill-Spill-Merge algorithm that
utilizes our depression hierarchy to rapidly process and distribute runoff.
Runoff fills depressions, which then overflow and spill into their neighbors. If
both a depression and its neighbor fill, they merge. We provide a detailed
explanation of the algorithm as well as results from two sample study areas. In
these case studies, the algorithm runs 90-2,600x faster (with a 2,000-63,000x
reduction in compute time) than the commonly-used Jacobi iteration and produces
a more accurate output. Complete, well-commented, open-source code is available
on Github and Zenodo.


Prerequisites
-------------

Ensure you have a working compiler.

The following compilers are known to work: GCC7.5.0, GCC8.4.0, GCC9.3.0

The following compilers are known to be too old: GCC5.4.0

Although GDAL is not required to use the library, it is needed to run the
example program.

Install the prerequisites

### Linux

    sudo apt install libgdal-dev cmake

### Mac

    brew install gdal libomp cmake llvm

(GDAL should not be necessary for installation, but is recommended.)



Additional setup
----------------

### Mac

*Here we assume that you are using the Anaconda Python distribution.*

These additional possible requirements to compile Fill-Spill-Merge on MacOSX stem from the fact that those users with the [Anaconda Python distribution now comes packaged with its own compilers](https://www.anaconda.com/blog/utilizing-the-new-compilers-in-anaconda-distribution-5). This issue may exist in Linux as well, though due to the extensive and well-integrated `apt` repositories, Linux users are less likely to install Anaconda Python.

1. Download MacOSX10.9.sdk from https://github.com/phracker/MacOSX-SDKs and place it in /opt/MacOSX10.9.sdk
2. Point the cmake compiler at this SDK:
```
export CONDA_BUILD_SYSROOT=/opt/MacOSX10.9.sdk
```

Compilation
-----------

Be sure to acquire submodules either upon initially obtaining the repository:

    git clone --recurse-submodules -j8 https://github.com/r-barnes/Barnes2020-FillSpillMerge

Or afterwards by using the following within the repository itself:

    git submodule update --init --recursive

Afterwards, compile:

### Linux

    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release -DUSE_GDAL=ON ..
    make -j 4 # Set to number of CPUs for a faster compilation (4 here)

### Mac

    mkdir build
    cd build
    # Be sure to repoint the versioning in the following as necessary:
    cmake -D CMAKE_C_COMPILER="/usr/local/Cellar/llvm/10.0.0_3/bin/clang" \
          -D CMAKE_CXX_COMPILER="/usr/local/Cellar/llvm/10.0.0_3/bin/clang++" \
          -DUSE_GDAL=ON
          -DCMAKE_BUILD_TYPE=Release ..
    make -j 4 # Set to number of CPUs for a faster compilation (4 here)

### Additional options

 * Use `-DCODE_COVERAGE=ON` to enable code coverage report generation.
   When doing a build this will run all the tests and generate a coverage
   report.



Running Example Program
-----------------------

A test program, `fsm.exe` is generated by make. This program simulates pouring a
given amount of water onto every cell on a landscape and determining where it
all ends up. Running the program shows its command-line arguments.



Running Tests
-------------

Compiling as above will produce an executable `fsm_unittests.exe`. Running this
will perform all the tests associated with the code. Randomized testing is used
to cover a wide range of possible inputs both for individual functions (unit
testing) as well as the entirety of FillSpillMerge in an end-to-end fashion.

Fill-Spill-Merge
===========================================

**Title of Manuscript**:
 * Computing water flow through complex landscapes, Part 3: Fill-Spill-Merge: Flow routing in depression hierarchies (doi: TODO)

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



Compilation
-----------

Be sure to acquire submodules either upon initially obtaining the repository:

    git clone --recurse-submodules -j8 https://github.com/r-barnes/Barnes2020-FillSpillMerge

Or afterwards by using the following within the repository itself:

    git submodule update --recursive

Afterwards, compile:

    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release ..
    make -j 4 #Set to number of CPUs for a faster compilation

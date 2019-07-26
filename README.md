Depression Hierarchies and Fill-Spill-Merge
===========================================

**Title of Manuscript**:
 * Computing water flow through complex landscapes, Part 2: Finding hierarchies in depressions and morphological segmentations (doi: [10.5194/esurf-2019-34](https://doi.org/10.5194/esurf-2019-34))

**Authors**: Richard Barnes, Kerry Callaghan, Andrew Wickert

**Corresponding Author**: Richard Barnes (richard.barnes@berkeley.edu)

**Code Repositories**
 * [Author's GitHub Repository](https://github.com/r-barnes/Barnes2019-DepressionHierarchy)

This repository contains a reference implementation of the algorithms presented
in the manuscript above, along with information on acquiring the various
datasets used, and code to perform correctness tests.



Abstract
--------

Depressions – inwardly-draining regions of digital elevation models – present
difficulties for terrain analysis and hydrological modeling. Analogous
"depressions" also arise in image processing and morphological segmentation
where they may represent noise, features of interest, or both. Here we provide a
new data structure – the depression hierarchy – that captures the full topologic
and topographic complexity of depressions in a region. We treat depressions as
networks, in a way that is analogous to surface-water flow paths, in which
individual sub-depressions merge together to form meta-depressions in a process
that continues until they begin to drain externally. The hierarchy can be used
to selectively fill or breach depressions, or to accelerate dynamic models of
hydrological flow. Complete, well-commented, open-source code and correctness
tests are available on Github and Zenodo.



Compilation
-----------

Run

    mkdir build
    cd build
    cmake ..
    make -j 4 #Set to number of CPUs for a faster compilation

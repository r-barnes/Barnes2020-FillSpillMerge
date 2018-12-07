#include "Array2D.hpp"
#include "dephier_kerry.hpp"
#include "DisjointDenseIntSet.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <queue>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>

int main(int argc, char **argv){
  if(argc!=4){
    std::cerr<<"Syntax: "<<argv[0]<<" <Input> <Output> <OutGraph>"<<std::endl;
    return -1;
  }

  const std::string in_name   = argv[1];
  const std::string out_name  = argv[2];
  const std::string out_graph = argv[3];

  Array2D<float> dem(in_name,"value");   //Recharge (Percipitation minus Evapotranspiration)

  //Initialize labels to indicate that none of the cells are part of depressions
  Array2D<label_t> label   (dem.width(),dem.height(),NO_DEP);

  //Initialize flow directions to indicate that none of the cells flow anywhere
  Array2D<flowdir_t> flowdirs(dem.width(),dem.height(),NO_FLOW);

  //Label the ocean cells. This is a precondition for using
  //`GetDepressionHierarchy()`.
  #pragma omp parallel for
  for(int i=0;i<label.size();i++)
    if(dem(i)==0)
      label(i) = OCEAN;

  //Label all the depressions and get the hierarchy connecting them.
  
  //TODO: 0. Calculate the number of cells within and volume of each depression.
  //This will take place inside of GetDepressionHierarchy.
  const auto deps = GetDepressionHierarchy<float,Topology::D8>(dem, label, flowdirs);

  //TODO 1. DONE. Get flow directions for all cells. 

  //TODO 2. Perform a flow accumulation moving water downhill and filling
  //groundwater as you go. Look at `misc/main.cpp`. Recall that we did this by
  //counting depressions, finding peaks, and using a queue to control a breadth-
  //first traversal from the peaks downwards.

  //TODO 3. As part of the above, when flow can't go downhill any farther, add
  //it to the `water_vol` for the appropriate depression. Use the labels array
  //to determine the appropriate depression.

  //TODO 4/5. Perform a depth-first post-order traversal of the depression
  //hierarchy (start with depressions for which `parent==NO_PARENT`. When you
  //reach the leaves if `water_vol>dep_vol` then try to overflow into the
  //neighbouring depression if its `water_vol<dep_vol`. To do so search
  //neighbour cells of `out_cell` for the lowest cell labeled `odep` and follow
  //that one's flow path until it terminates. Add excess water to that
  //depression's `water_vol`. After both child depressions have been visited
  //they will be finished trying to share their water and their excess water is
  //added to their parent's `water_vol`. Repeat the overflow attempt.

  //TODO 6: Adjust the hydrologic elevations. If a depression has `water_vol>0`
  //then all of its child depression are full. What remains is to use a
  //priority-queue and the Water Level Equation (see dephier.cpp) to find which
  //cells should be flooded. If a depression is a leaf depression (no children)
  //then start the PQ at the pit_cell. Otherwise, start at the pit cell of any
  //child depression since all children are flooded to at least the level of the
  //`out_elev` connecting the uppermost two children in the hierarchy. Cells
  //should be added to the queue with elevation equal to at least the `out_elev`
  //of the depression's children and, if `water_vol<dep_vol` should not exceed
  //the `out_elev` of the depression itself.



  //TODO: Remove. For viewing test cases.
  if(label.width()<1000){
    for(int y=0;y<label.height();y++){
      for(int x=0;x<label.width();x++)
        std::cout<<std::setw(3)<<label(x,y)<<" ";
      std::cout<<std::endl;
    }

    //GraphViz dot-style output for drawing depression hierarchy graphs.
    std::ofstream fgraph(out_graph);
    fgraph<<"digraph {\n";
    for(unsigned int i=0;i<deps.size();i++)
      fgraph<<i<<" -> "<<deps[i].parent<<";\n";
    fgraph<<"}\n";
  }

  SaveAsNetCDF(dem,out_name+"-dem.nc","value");
  SaveAsNetCDF(label,out_name+"-labels_raw.nc","value");

  LastLayer(label, dem, deps);

  SaveAsNetCDF(label,out_name+"-labels_proc.nc","value");

  return 0;
}


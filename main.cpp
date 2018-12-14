#include "debugging_utilities.hpp"
#include "dh_flow.hpp"
#include "../common/netcdf.hpp"
#include <iostream>
#include <richdem/common/Array2D.hpp>
#include <string>


int main(int argc, char **argv){
  const float  OCEAN_LEVEL = 0;  //ocean_level in the topo file must be lower than any non-ocean cell. 
  
  if(argc!=4){
    std::cerr<<"Syntax: "<<argv[0]<<" <Input> <Output> <OutGraph>"<<std::endl;
    return -1;
  }

  rd::Timer timer_overall;
  rd::Timer timer_io;
  timer_overall.start();

  const std::string in_name   = argv[1];
  const std::string out_name  = argv[2];
  const std::string out_graph = argv[3];

  timer_io.start();
  rd::Array2D<float> topo = LoadData<float>(in_name,std::string("value"));   //Recharge (Percipitation minus Evapotranspiration)
  timer_io.stop();

  PrintDEM("Topography", topo);

  rd::Array2D<float>     wtd     (topo.width(), topo.height(), 1      ); //All cells have some water
  rd::Array2D<label_t>   label   (topo.width(), topo.height(), NO_DEP ); //No cells are part of a depression
  rd::Array2D<flowdir_t> flowdirs(topo.width(), topo.height(), NO_FLOW); //No cells flow anywhere

  //Label the ocean cells. This is a precondition for using
  //`GetDepressionHierarchy()`.
  #pragma omp parallel for
  for(unsigned int i=0;i<label.size();i++)
    if(topo.isNoData(i) || topo(i)==OCEAN_LEVEL){ //Ocean Level is assumed to be lower than any other cells (even Death Valley)
      label(i) = OCEAN;
      wtd  (i) = 0;
    }

  //Generate flow directions, label all the depressions, and get the hierarchy
  //connecting them
  auto deps = GetDepressionHierarchy<float,Topology::D8>(topo, label, flowdirs);

  PrintDEM("labels", label);

  PrintDepressionInfo(deps);

  //TODO: Remove. For viewing test cases.
  if(label.width()<1000){
    //GraphViz dot-style output for drawing depression hierarchy graphs.
    std::ofstream fgraph(out_graph);
    fgraph<<"digraph {\n";
    for(int i=0;i<(int)deps.size();i++){
      fgraph<<i<<" -> "<<deps[i].parent;
      if(deps[i].parent!=NO_VALUE && (deps[i].parent==OCEAN || !(deps[deps[i].parent].lchild==i || deps[deps[i].parent].rchild==i)))
        fgraph<<" [color=\"blue\"]";
      fgraph<<";\n";
    }
    fgraph<<"}\n";
  }

  SaveAsNetCDF(topo, out_name+"-topo.nc",       "value");
  SaveAsNetCDF(label,out_name+"-labels_raw.nc", "value");
  SaveAsNetCDF(label,out_name+"-labels_proc.nc","value");




  SurfaceWater(topo, wtd, label,deps,flowdirs);


  PrintDEM("wtd", wtd);




  std::unordered_map<label_t, label_t> jump_table;
  Overflow(OCEAN, deps, jump_table);
  jump_table = std::unordered_map<label_t, label_t>();

  PrintDepressionInfo(deps);

  //Sanity checks
  for(int d=1;d<(int)deps.size();d++){
    const auto &dep = deps.at(d);
    assert(dep.water_vol==0 || dep.water_vol<=dep.dep_vol);
    assert(dep.water_vol==0 || (dep.lchild==NO_VALUE && dep.rchild==NO_VALUE) || (dep.lchild!=NO_VALUE && deps.at(dep.lchild).water_vol<dep.water_vol));
    assert(dep.water_vol==0 || (dep.lchild==NO_VALUE && dep.rchild==NO_VALUE) || (dep.rchild!=NO_VALUE && deps.at(dep.rchild).water_vol<dep.water_vol));
  }

  PrintDEM("wtd", wtd, 9);

  // for(auto &depression:deps){//(unsigned int d=0;d<deps.size();d++){
  //   std::cerr<<"Here's the list of all depressions with their parents and children: "
  //            <<depression.dep_label<<" "
  //            <<std::setw(3)<<depression.parent   <<" "
  //            <<std::setw(3)<<depression.lchild   <<" "
  //            <<std::setw(3)<<depression.rchild   <<" "
  //            <<depression.water_vol
  //            <<std::endl;
  //   std::cerr<<"\tOcean-linked = ";
  //   for(auto x: depression.ocean_linked)
  //     std::cerr<<x<<" ";
  //   std::cerr<<std::endl;
  // }
  

  std::cerr<<"\n\n\033[91m#######################Finding Filled\033[39m"<<std::endl;
  Find_filled(OCEAN,deps,topo,label,wtd);                              //This should check everything that is an immediate child of the ocean, so we're supposed to hit all the depressions like this. 
       
  SaveAsNetCDF(wtd,out_name+"-wtd.nc","value");

  for(int i=0;i<topo.size();i++)
    if(!topo.isNoData(i))
      wtd(i) += topo(i);

  SaveAsNetCDF(wtd,out_name+"-flooded.nc","value");

  rd::FillDepressions<rd::Topology::D8>(topo);
  SaveAsNetCDF(topo,out_name+"-filled.nc","value");



  std::cerr<<"Finished"<<std::endl;

  if(topo.width()<1000){
    PrintDEM("topo",     topo    );
    PrintDEM("Flowdirs", flowdirs);
    PrintDEM("wtd",      wtd     );
    PrintDEM("labels",   label   );
  }

  std::cerr<<"Wall-time = "<<timer_overall.stop()  <<" s"<<std::endl;
  std::cerr<<"IO time   = "<<timer_io.accumulated()<<" s"<<std::endl;

  return 0;
}

#include "dephier.hpp"

#include <richdem/common/Array2D.hpp>
#include <richdem/depressions/Barnes2014.hpp>
#include <richdem/depressions/Zhou2016.hpp>
#include <richdem/depressions/Wei2018.hpp>

#include <iostream>
#include <string>
#include <stdexcept>

namespace rd = richdem;
namespace dh = richdem::dephier;

int main(int argc, char **argv){
  if(argc!=4){
    std::cout<<"Syntax: "<<argv[0]<<" <Input> <Ocean Level> <Output Prefix>"<<std::endl;
    return -1;
  }

  const std::string in_name     = argv[1];
  const float       ocean_level = std::stod(argv[2]);
  const std::string out_prefix = argv[3];

  std::cout<<"m Input file  = "<<argv[1]<<std::endl;
  std::cout<<"m Ocean level = "<<argv[2]<<std::endl;
  std::cout<<"m Out prefix  = "<<argv[3]<<std::endl;

  rd::Timer timer_io;
  timer_io.start();
  const rd::Array2D<float> topo_orig(in_name);
  timer_io.stop();
  std::cout<<"IO time = "<<timer_io.accumulated()<<" s"<<std::endl;

  std::cout<<"m Data width  = "<<topo_orig.width ()<<std::endl;
  std::cout<<"m Data height = "<<topo_orig.height()<<std::endl;
  std::cout<<"m Data cells  = "<<topo_orig.numDataCells()<<std::endl;

  dh::DepressionHierarchy<float> deps;
  rd::Array2D<dh::dh_label_t> label   (topo_orig.width(), topo_orig.height(), dh::NO_DEP ); //No cells are part of a depression
  rd::Array2D<rd::flowdir_t>  flowdirs(topo_orig.width(), topo_orig.height(), rd::NO_FLOW); //No cells flow anywhere

  {
    std::cout<<"p Labeling the ocean cells..."<<std::endl;
    #pragma omp parallel for
    for(unsigned int i=0;i<label.size();i++){
      if(topo_orig.isNoData(i) || topo_orig(i)==ocean_level){ //Ocean Level is assumed to be lower than any other cells (even Death Valley)
        label(i) = dh::OCEAN;
      }
    }
  }

  {
    rd::Timer timer;
    auto topo2 = topo_orig;
    timer.start();
    rd::PriorityFlood_Barnes2014<rd::Topology::D8>(topo2);
    timer.stop();
    std::cout<<"Barnes2014 depresison filling time = "<<timer.accumulated()<<std::endl;
  }

  {
    rd::Timer timer;
    auto topo2 = topo_orig;
    timer.start();
    rd::PriorityFlood_Zhou2016(topo2);
    timer.stop();
    std::cout<<"Zhou2016 depresison filling time = "<<timer.accumulated()<<std::endl;
  }

  {
    rd::Timer timer;
    auto topo2 = topo_orig;
    timer.start();
    rd::PriorityFlood_Wei2018(topo2);
    timer.stop();
    std::cout<<"Wei2018 depresison filling time = "<<timer.accumulated()<<std::endl;
  }

  {
    rd::Timer timer;
    auto topo2 = topo_orig;
    timer.start();
    deps = dh::GetDepressionHierarchy<float,rd::Topology::D8>(topo2, label, flowdirs);
    timer.stop();
    std::cout<<"Depression Hierarchy time = "<<timer.accumulated()<<std::endl;
  }

  //TODO: Remove. For viewing test cases.
  // {
  //   //GraphViz dot-style output for drawing depression hierarchy graphs.
  //   std::ofstream fgraph(out_graph);
  //   fgraph<<"digraph {\n";
  //   for(unsigned int i=0;i<deps.size();i++){
  //     fgraph<<i<<" -> "<<deps[i].parent;
  //     if(deps[i].parent!=dh::NO_VALUE && (deps[i].parent==dh::OCEAN || !(deps[deps[i].parent].lchild==i || deps[deps[i].parent].rchild==i)))
  //       fgraph<<" [color=\"blue\"]";
  //     fgraph<<";\n";
  //   }
  //   fgraph<<"}\n";
  // }

  {
    std::cout<<"p Generating depression stats..."<<std::endl;

    std::ofstream fout(out_prefix+"-depstats.csv");
    fout<<"pit_elev,out_elev,cell_count,dep_vol,water_vol,marginal_vol\n";
    for(const auto &dep: deps){
      auto marginal_vol = dep.dep_vol;
      if(dep.lchild!=dh::NO_VALUE)
        marginal_vol -= deps.at(dep.lchild).dep_vol;
      if(dep.rchild!=dh::NO_VALUE)
        marginal_vol -= deps.at(dep.rchild).dep_vol;
      fout<<dep.pit_elev
          <<" "<<dep.out_elev
          <<" "<<dep.cell_count
          <<" "<<dep.dep_vol
          <<" "<<dep.water_vol
          <<" "<<marginal_vol
          <<"\n";
    }
  }

  {
    std::cout<<"p Saving base labels..."<<std::endl;

    label.saveGDAL(out_prefix+"-baselabels.tif");
  }

  {
    std::cout<<"p Generating watersheds..."<<std::endl;

    auto label2 = label;
    #pragma omp parallel for collapse(2)
    for(int y=0;y<label2.height();y++)
    for(int x=0;x<label2.width();x++){
      auto mylabel = label2(x,y);
      while(mylabel!=dh::OCEAN && mylabel!=dh::NO_PARENT && deps.at(mylabel).parent!=dh::OCEAN && deps.at(mylabel).parent!=dh::NO_PARENT){
        mylabel = deps.at(mylabel).parent;
      }
      label2(x,y) = mylabel;
    }

    label2.saveGDAL(out_prefix+"-toplabels.tif");
  }

  //Filter depressions by relabeling small ones with their parent's data
  {
    std::cout<<"p Filtering small depressions..."<<std::endl;

    auto label2 = label;
    #pragma omp parallel for collapse(2)
    for(int y=0;y<label2.height();y++)
    for(int x=0;x<label2.width();x++){
      auto mylabel = label2(x,y);
      while(
           mylabel!=dh::OCEAN 
        && deps.at(mylabel).parent!=dh::OCEAN
        && deps.at(mylabel).cell_count<30
      ){
        mylabel = deps.at(mylabel).parent;
      }
      label2(x,y) = mylabel;
    }

    label2.saveGDAL(out_prefix+"-filteredlabels.tif");
  }


  return 0;
}

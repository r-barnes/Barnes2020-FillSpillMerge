#include "fill_spill_merge.hpp"
#include "netcdf.hpp"
#include <iostream>
#include <richdem/common/Array2D.hpp>
#include "priority_flood.hpp"
#include <string>
#include <stdexcept>

namespace rd = richdem;
namespace dh = richdem::dephier;

int main(int argc, char **argv){
  if(argc!=4){
    std::cout<<"Syntax: "<<argv[0]<<" <Input> <Out prefix> <Ocean Level>"<<std::endl;
    return -1;
  }

  const std::string in_name     = argv[1];
  const std::string out_name    = argv[2];
  const float       ocean_level = std::stod(argv[3]);

  std::cout<<"m Processing    = "<<argv[1]<<std::endl;
  std::cout<<"m Output prefix = "<<argv[2]<<std::endl;
  std::cout<<"m Ocean level   = "<<argv[3]<<std::endl;

  rd::Timer timer_io;
  timer_io.start();
  rd::Array2D<float> topo(in_name);   //Recharge (Percipitation minus Evapotranspiration)
  timer_io.stop();

  std::cout<<"m Data width  = "<<topo.width ()<<std::endl;
  std::cout<<"m Data height = "<<topo.height()<<std::endl;
  std::cout<<"m Data cells  = "<<topo.numDataCells()<<std::endl;

  rd::Array2D<float>          wtd     (topo.width(), topo.height(), 10000       ); //All cells have some water
  rd::Array2D<dh::dh_label_t> label   (topo.width(), topo.height(), dh::NO_DEP );  //No cells are part of a depression
  rd::Array2D<rd::flowdir_t>  flowdirs(topo.width(), topo.height(), rd::NO_FLOW);  //No cells flow anywhere

  wtd.setNoData(topo.noData());

  //Label the ocean cells. This is a precondition for using
  //`GetDepressionHierarchy()`.
  #pragma omp parallel for
  for(unsigned int i=0;i<label.size();i++){
    if(topo.isNoData(i) || topo(i)==ocean_level){ //Ocean Level is assumed to be lower than any other cells (even Death Valley)
      label(i) = dh::OCEAN;
      wtd  (i) = 0;
    }
  }

  //Generate flow directions, label all the depressions, and get the hierarchy
  //connecting them
  rd::Timer timer_dh;
  timer_dh.start();
  auto deps = dh::GetDepressionHierarchy<float,rd::Topology::D8>(topo, label, flowdirs);
  timer_dh.stop();
  std::cout<<"t DH time  = "<<timer_dh.accumulated()<<" s"<<std::endl;

  rd::Timer timer_fsm;
  timer_fsm.start();
  dh::FillSpillMerge(topo, label, flowdirs, deps, wtd, ocean_level);
  timer_fsm.stop();
  std::cout<<"t FSM time = "<<timer_fsm.accumulated()<<" s"<<std::endl;

  for(unsigned int i=0;i<topo.size();i++)
    if(!topo.isNoData(i))
      wtd(i) += topo(i);

  SaveAsNetCDF(wtd,out_name+"-flooded.nc","value");
  rd::PriorityFlood_Barnes2014_OceanInit<rd::Topology::D8>(topo, ocean_level);
  SaveAsNetCDF(topo,out_name+"-filled.nc","value");

  double max_diff = 0;
  unsigned count  = 0;
  double avg_diff = 0;
  for(unsigned int i=0;i<topo.size();i++){
    if(topo.isNoData(i))
      continue;
    const auto diff = static_cast<double>(std::abs(wtd(i)-topo(i)));
    max_diff=std::max(diff*diff,max_diff);
    count++;
    avg_diff+=(diff-avg_diff)/count;
  }
  // SaveAsNetCDF(diff,out_name+"-diff.nc","value");

  std::cout<<"m Max diff = "<<std::sqrt(max_diff)<<std::endl;
  std::cout<<"m Avg diff = "<<avg_diff<<std::endl;

  std::cout<<"p Finished"<<std::endl;
  std::cout<<"t IO time  = "<<timer_io.accumulated()<<" s"<<std::endl;

  return 0;
}

#include <fsm/fill_spill_merge.hpp>
#include "netcdf.hpp"
#include <iostream>
#include <richdem/common/Array2D.hpp>
#include <string>
#include <stdexcept>

namespace rd = richdem;
namespace dh = richdem::dephier;

int main(int argc, char **argv){
  if(argc!=5){
    std::cout<<"Syntax: "<<argv[0]<<" <Input> <Runoff> <Ocean Level> <Outfile>"<<std::endl;
    return -1;
  }

  const std::string in_name     = argv[1];
  const double      flood_level = std::stod(argv[2]);
  const float       ocean_level = std::stod(argv[3]);
  const std::string out_name    = argv[4];

  std::cout<<"m Processing  = "<<in_name    <<std::endl;
  std::cout<<"m Runoff      = "<<flood_level<<std::endl;
  std::cout<<"m Ocean level = "<<ocean_level<<std::endl;
  std::cout<<"m Out name    = "<<out_name   <<std::endl;

  rd::Timer timer_overall;
  timer_overall.start();
  rd::Array2D<float> topo(in_name);   //Recharge (Percipitation minus Evapotranspiration)

  rd::Array2D<float>          wtd     (topo.width(), topo.height(), flood_level); //All cells have some water
  rd::Array2D<dh::dh_label_t> label   (topo.width(), topo.height(), dh::NO_DEP );  //No cells are part of a depression
  rd::Array2D<rd::flowdir_t>  flowdirs(topo.width(), topo.height(), rd::NO_FLOW);  //No cells flow anywhere

  wtd.setNoData(topo.noData());

  //Label the ocean cells. This is a precondition for using
  //`GetDepressionHierarchy()`.

  int ocean_level_exists = 0;
  #pragma omp parallel for
  for(unsigned int i=0;i<label.size();i++){
    if(topo.isNoData(i) || topo(i)==ocean_level){ //Ocean Level is assumed to be lower than any other cells (even Death Valley)
      label(i) = dh::OCEAN;
      wtd  (i) = 0;
    }
    if(topo(i) == ocean_level)
      ocean_level_exists = 1;
  }

  if(ocean_level_exists == 0)
    throw std::runtime_error("There are no ocean_level cells in this topography!");

  //Generate flow directions, label all the depressions, and get the hierarchy
  //connecting them
  auto deps = dh::GetDepressionHierarchy<float,rd::Topology::D8>(topo, label, flowdirs);
  dh::FillSpillMerge(topo, label, flowdirs, deps, wtd, ocean_level);

  for(unsigned int i=0;i<topo.size();i++)
    if(!topo.isNoData(i))
      wtd(i) += topo(i);

  SaveAsNetCDF(wtd,out_name+"-flooded.nc","value");

  timer_overall.stop();
  std::cout<<"t Total time  = "<<timer_overall.accumulated()<<" s"<<std::endl;

  return 0;
}

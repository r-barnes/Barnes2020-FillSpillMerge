#include <fsm/fill_spill_merge.hpp>
#include <richdem/common/Array2D.hpp>
#include <richdem/misc/misc_methods.hpp>
#include <richdem/common/gdal.hpp>

#include <iostream>
#include <string>
#include <stdexcept>

namespace rd = richdem;
namespace dh = richdem::dephier;

int main(int argc, char **argv){
  if(argc!=5){
    std::cout<<"Pours a given amount of water onto every cell of a landscape and determines where it all ends up\n"<<std::endl;
    std::cout<<"Syntax: "<<argv[0]<<" <Input> <Output> <Surface Water Level> <Ocean Level>"<<std::endl;
    std::cout<<"Ocean cells are detected at the DEM perimeter and added from there"<<std::endl;
    return -1;
  }

  const std::string in_name          = argv[1];
  const std::string out_name         = argv[2];
  const double      surf_water_level = std::stod(argv[3]);
  const double      ocean_level      = std::stod(argv[4]);

  std::cout<<"m Input DEM           = "<<argv[1]<<std::endl;
  std::cout<<"m Output prefix       = "<<argv[2]<<std::endl;
  std::cout<<"m Surface water level = "<<argv[3]<<std::endl;
  std::cout<<"m Ocean level         = "<<argv[4]<<std::endl;

  rd::Timer timer_io;
  timer_io.start();
  rd::Array2D<double> topo(in_name);
  timer_io.stop();

  std::cout<<"m Data width  = "<<topo.width ()<<std::endl;
  std::cout<<"m Data height = "<<topo.height()<<std::endl;
  std::cout<<"m Data cells  = "<<topo.numDataCells()<<std::endl;

  rd::Array2D<double>          wtd     (topo.width(), topo.height(), surf_water_level); //All cells have some water
  rd::Array2D<dh::dh_label_t> label   (topo.width(), topo.height(), dh::NO_DEP );      //No cells are part of a depression
  rd::Array2D<rd::flowdir_t>  flowdirs(topo.width(), topo.height(), rd::NO_FLOW);      //No cells flow anywhere

  wtd.setNoData(topo.noData());

  //Label the ocean cells. This is a precondition for using
  //`GetDepressionHierarchy()`.
  rd::BucketFillFromEdges<rd::Topology::D8>(topo, label, ocean_level, dh::OCEAN);

  //Make NoData cells also ocean cells. Ocean has no water on it to begin with.
  #pragma omp parallel for
  for(unsigned int i=0;i<label.size();i++){
    if(topo.isNoData(i) || label(i)==dh::OCEAN){
      label(i) = dh::OCEAN;
      wtd  (i) = 0;
    }
  }

  rd::Timer timer_calc;
  timer_calc.start();

  //Generate flow directions, label all the depressions, and get the hierarchy
  //connecting them
  auto deps = dh::GetDepressionHierarchy<double,rd::Topology::D8>(topo, label, flowdirs);

  dh::FillSpillMerge(topo, label, flowdirs, deps, wtd);

  timer_calc.stop();

  timer_io.start();

  //Output the water table depth
  wtd.saveGDAL(out_name+"-wtd.nc");

  for(unsigned int i=0;i<topo.size();i++)
    if(!topo.isNoData(i))
      wtd(i) += topo(i);

  //Output the new height of the hydraulic surface
  wtd.saveGDAL(out_name+"-hydraulic-surface-height.nc");

  timer_io.stop();

  std::cout<<"Finished."<<std::endl;
  std::cout<<"IO time   = "<<timer_io.accumulated()  <<" s"<<std::endl;
  std::cout<<"Calc time = "<<timer_calc.accumulated()<<" s"<<std::endl;

  return 0;
}

#include <dephier/dephier.hpp>
#include <fsm/fill_spill_merge.hpp>
#include <richdem/common/loaders.hpp>
#include <richdem/terrain_generation.hpp>

#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

using namespace richdem;
using namespace richdem::dephier;

const double ocean_level = -1;



void time_fsm(std::mt19937_64 &gen, const int size){
  static std::uniform_int_distribution<uint32_t> seed_dist;

  auto dem = perlin(size, seed_dist(gen));

  Array2D<double>     wtd     (dem.width(),dem.height(), 10);
  Array2D<dh_label_t> labels  (dem.width(),dem.height(), NO_DEP);
  Array2D<flowdir_t>  flowdirs(dem.width(),dem.height(), NO_FLOW);

  dem.setEdges(ocean_level);
  labels.setEdges(OCEAN);

  Timer time_dh, time_fsm, time_both;

  time_both.start();
  time_dh.start();
  auto deps = GetDepressionHierarchy<double,Topology::D8>(dem, labels, flowdirs);
  time_dh.stop();

  time_fsm.start();
  FillSpillMerge<double,double>(dem, labels, flowdirs, deps, wtd);
  time_fsm.stop();
  time_both.stop();

  std::cout<<"DH,  "<<std::setw(5)<<size<<","<<time_dh.accumulated()<<std::endl;
  std::cout<<"FSM, "<<std::setw(5)<<size<<","<<time_fsm.accumulated()<<std::endl;
  std::cout<<"Both,"<<std::setw(5)<<size<<","<<time_both.accumulated()<<std::endl;
}



int main(){
  std::mt19937_64 gen;

  //Sizes to test at
  const std::vector<int> sizes = {10,50,100,500,1000,5000,10000,50000};

  //Number of times to run each test
  const int reps = 3;

  for(const auto size: sizes)
  for(int i=0;i<reps;i++)
    time_fsm(gen, size);

  return 0;
}
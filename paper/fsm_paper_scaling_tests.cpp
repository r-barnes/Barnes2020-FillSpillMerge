#include <richdem/terrain_generation.hpp>
#include <richdem/common/loaders.hpp>
#include <fsm/fill_spill_merge.hpp>
#include <dephier/dephier.hpp>

namespace rd = richdem;
namespace dh = richdem::dephier;

/**
  @brief Pads a 2D array on all sides by cells of a given value

  @param[in] pad_width Positive number of cells in the padding
  @param[in] value     Value to give to the padding cells

  @return A new array based on the old one
*/
rd::Array2D<double> pad(const rd::Array2D<double> &inp, const int pad_width, const double value){
  rd::Array2D<double> temp(inp.width()+2*pad_width, inp.height()+2*pad_width);

  temp.setAll(value);

  if(pad_width<0)
    throw std::runtime_error("Cannot have a negative padding!");

  for(size_t y=0;y<inp.height();y++)
  for(size_t x=0;x<inp.width();x++){
    temp(x+pad_width,y+pad_width) = inp(x,y);
  }

  return temp;
}



int main(){
  auto dem = richdem::perlin(2000);

  const double ocean_level = -1;

  dem = pad(dem, 1, ocean_level);

  //Add water everywhere. Since elevation of the random terrain is in the range
  //[0,1], this is actually quite a bit.

  rd::Array2D<double>         wtd     (dem.width(),dem.height(), 10);
  rd::Array2D<dh::dh_label_t> labels  (dem.width(),dem.height(), dh::NO_DEP);
  rd::Array2D<rd::flowdir_t>  flowdirs(dem.width(),dem.height(), rd::NO_FLOW);

  labels.setAll(dh::NO_DEP);
  for(int y=0;y<labels.height();y++)
  for(int x=0;x<labels.width() ;x++){
    if(labels.isEdgeCell(x,y)){
      labels(x,y) = dh::OCEAN;
      wtd(x,y)    = 0;
    }
  }

  flowdirs.setAll(rd::NO_FLOW);

  auto deps = dh::GetDepressionHierarchy<double,rd::Topology::D8>(dem, labels, flowdirs);
  dh::FillSpillMerge<double,double>(dem, labels, flowdirs, deps, wtd, ocean_level);

  SaveGDAL(dem, "/z/out.tif");
}
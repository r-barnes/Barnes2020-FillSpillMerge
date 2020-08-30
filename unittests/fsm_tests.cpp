#include "doctest.h"
#include <fsm/fill_spill_merge.hpp>

#include <richdem/terrain_generation.hpp>

#include <random>

using namespace richdem;
using namespace richdem::dephier;

template<class T>
double MaxArrayDiff(const Array2D<T> &a, const Array2D<T> &b){
  double max_diff = 0;
  for(auto i=a.i0();i<a.size();i++){
    max_diff = std::max(max_diff, (double)std::abs(a(i)-b(i)));
  }
  return max_diff;
}

template<class T>
bool ArrayValuesEqual(const Array2D<T> &a, const Array2D<T> &b){
  for(auto i=a.i0();i<a.size();i++){
    if(a(i)!=b(i))
      return false;
  }
  return true;
}

template<class T>
bool ArrayValuesAllEqual(const Array2D<T> &a, const T val){
  for(auto i=a.i0();i<a.size();i++){
    if(a(i)!=val)
      return false;
  }
  return true;
}


TEST_CASE("Depression volume"){
  CHECK(DepressionVolume(2, 5, 10)==0);
  CHECK(DepressionVolume(3, 5, 10)==5);
  CHECK(DepressionVolume(4, 5, 10)==10);
}

TEST_CASE("Determine water level"){
  SUBCASE("Depression volume exactly equals water volume"){
    double sill_wtd = -2;
    const auto water_level = DetermineWaterLevel(sill_wtd, 10, 4, 5, 10);
    CHECK(sill_wtd==-2);
    CHECK(water_level==4); //Water elevation equals the sill elevation
  }

  SUBCASE("Water volume is less than the depression volume"){
    double sill_wtd = -2;
    const auto water_level = DetermineWaterLevel(sill_wtd, 8, 4, 5, 10);
    CHECK(sill_wtd==-2);
    CHECK(water_level==18/5.0); //Water elevation equals the sill elevation
  }

  SUBCASE("Water volume is greater than the depression volume"){
    double sill_wtd = -2;
    const auto water_level = DetermineWaterLevel(sill_wtd, 12, 4, 5, 10);
    CHECK(sill_wtd==0);
    //Water elevation equals the sill elevation since the sill absorbs all
    //excess
    CHECK(water_level==4);
  }
}



TEST_CASE("MoveWaterIntoPits 1"){
  const Array2D<double> topo = {
      {-9, -9, -9, -9, -9, -9, -9, -9, -9, -9},
      {-9,  9,  9,  9,  9,  9,  9,  9,  9, -9},
      {-9,  9,  8,  8,  8,  8,  7,  6,  9, -9},
      {-9,  9,  8,  7,  8,  7,  6,  5,  9, -9},
      {-9,  9,  8,  7,  8,  6,  5,  4,  9, -9},
      {-9,  9,  8,  8,  8,  5,  4,  3,  9, -9},
      {-9,  9,  7,  6,  5,  4,  3,  2,  9, -9},
      {-9,  9,  7,  6,  5,  4,  3,  1,  9, -9},
      {-9,  9,  9,  9,  9,  9,  9,  9,  9, -9},
      {-9, -9, -9, -9, -9, -9, -9, -9, -9, -9},
  };

  Array2D<dh_label_t> label(topo.width(), topo.height(), NO_DEP);
  for(int y=0;y<topo.height();y++)
  for(int x=0;x<topo.width();x++){
    if(topo.isEdgeCell(x,y))
      label(x,y) = OCEAN;
  }

  Array2D<flowdir_t> flowdirs(topo.width(), topo.height(), NO_FLOW);
  Array2D<double> wtd(topo.width(), topo.height(), 0);

  auto DH = GetDepressionHierarchy<double,Topology::D8>(topo, label, flowdirs);

  const Array2D<dh_label_t> label_good = {
    {0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0},
    {0,0,2,2,1,1,1,1,0,0},
    {0,0,2,2,1,1,1,1,0,0},
    {0,0,2,2,1,1,1,1,0,0},
    {0,0,1,1,1,1,1,1,0,0},
    {0,0,1,1,1,1,1,1,0,0},
    {0,0,1,1,1,1,1,1,0,0},
    {0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0}
  };

  const Array2D<flowdir_t> flowdirs_good = {
    {0,0,0,0,0,0,0,0,0,0},
    {0,8,4,4,4,4,4,4,6,0},
    {0,8,6,7,6,6,6,7,6,0},
    {0,8,6,7,6,6,6,7,6,0},
    {0,8,5,0,6,6,6,7,6,0},
    {0,8,6,6,6,6,6,7,6,0},
    {0,8,6,5,6,5,6,7,6,0},
    {0,8,5,4,5,4,5,0,6,0},
    {0,6,6,6,6,6,6,6,6,0},
    {0,0,0,0,0,0,0,0,0,0}
  };

  CHECK(ArrayValuesEqual(label,label_good));
  CHECK(ArrayValuesEqual(flowdirs,flowdirs_good));

  wtd.setAll(1);

  MoveWaterIntoPits<double, double>(topo, label, flowdirs, DH, wtd);

  CHECK(ArrayValuesAllEqual(wtd,0.0));

  CHECK(DH.at(0).water_vol==64);
  CHECK(DH.at(1).water_vol==30);
  CHECK(DH.at(2).water_vol== 6);

  CHECK(DH.at(0).parent==NO_VALUE);
  CHECK(DH.at(1).parent==3);
  CHECK(DH.at(2).parent==3);
  CHECK(DH.at(3).parent==0);

  CHECK(std::isnan(DH.at(0).dep_vol));
  CHECK(DH.at(1).dep_vol== 73);
  CHECK(DH.at(2).dep_vol==  2);
  CHECK(DH.at(3).dep_vol==111);
}
//TODO: Add second test case with more tests and clearer outlets



TEST_CASE("Backfill Depression"){
  const Array2D<double> topo = {
      {-9, -9, -9, -9, -9, -9, -9, -9, -9, -9},
      {-9,  6,  6,  6,  6,  6,  6,  6,  6, -9},
      {-9,  6,  1,  6,  1,  1,  6,  1,  6, -9},
      {-9,  6,  1,  6,  1,  3,  6,  1,  6, -9},
      {-9,  6,  1,  6,  2,  1,  4,  1,  6, -9},
      {-9,  6,  1,  6,  1,  1,  6,  1,  6, -9},
      {-9,  6,  1,  6,  6,  6,  6,  1,  6, -9},
      {-9,  6,  1,  1,  1,  1,  1,  1,  6, -9},
      {-9,  6,  6,  6,  6,  6,  6,  6,  6, -9},
      {-9, -9, -9, -9, -9, -9, -9, -9, -9, -9},
  };

  Array2D<double> wtd(topo.width(), topo.height());

  std::vector<flat_c_idx> cells_affected = {
    topo.xyToI(4,2), topo.xyToI(5,2),
    topo.xyToI(4,3), topo.xyToI(5,3),
    topo.xyToI(4,4), topo.xyToI(5,4),
    topo.xyToI(4,5), topo.xyToI(5,5),
  };

  BackfillDepression(4.0, topo, wtd, cells_affected);

  const Array2D<double> wtd_good = {
      { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
      { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
      { 0,  0,  0,  0,  3,  3,  0,  0,  0,  0},
      { 0,  0,  0,  0,  3,  1,  0,  0,  0,  0},
      { 0,  0,  0,  0,  2,  3,  0,  0,  0,  0},
      { 0,  0,  0,  0,  3,  3,  0,  0,  0,  0},
      { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
      { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
      { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
      { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
  };

  CHECK(ArrayValuesEqual(wtd,wtd_good));
}

TEST_CASE("FillDepressions"){
  const Array2D<double> topo = {
      {-9, -9, -9, -9, -9, -9, -9, -9, -9, -9},
      {-9,  6,  6,  6,  6,  6,  6,  6,  6, -9},
      {-9,  6,  1,  6,  1,  1,  6,  1,  6, -9},
      {-9,  6,  1,  6,  3,  3,  6,  1,  6, -9},
      {-9,  6,  1,  6,  2,  1,  4,  1,  6, -9},
      {-9,  6,  1,  6,  1,  1,  6,  1,  6, -9},
      {-9,  6,  1,  6,  6,  6,  6,  1,  6, -9},
      {-9,  6,  1,  1,  1,  1,  1,  1,  6, -9},
      {-9,  6,  6,  6,  6,  6,  6,  6,  6, -9},
      {-9, -9, -9, -9, -9, -9, -9, -9, -9, -9},
  };

  const Array2D<dh_label_t> label = {
      { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
      { 0,  1,  1,  2,  2,  2,  2,  1,  1,  0},
      { 0,  1,  1,  2,  2,  2,  2,  1,  1,  0},
      { 0,  1,  1,  2,  2,  2,  2,  1,  1,  0},
      { 0,  1,  1,  2,  2,  2,  2,  1,  1,  0},
      { 0,  1,  1,  2,  2,  2,  2,  1,  1,  0},
      { 0,  1,  1,  2,  2,  2,  2,  1,  1,  0},
      { 0,  1,  1,  1,  1,  1,  1,  1,  1,  0},
      { 0,  1,  1,  1,  1,  1,  1,  1,  1,  0},
      { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
  };

  Array2D<double> wtd(topo.width(), topo.height(), 0.0);

  const auto pit_cell = topo.xyToI(4,2);

  REQUIRE(topo(pit_cell)==1);

  const std::unordered_set<dh_label_t> dep_labels = {2};

  const double ocean_level = -9;

  SUBCASE("No water to add"){
    wtd(4,3) = -0.5;
    const auto wtd_good = wtd;
    const double water_vol = 0;
    FillDepressions(pit_cell, dep_labels, water_vol, topo, label, wtd, ocean_level);
    CHECK(wtd==wtd_good);
  }

  SUBCASE("Standard Case"){
    wtd.setAll(0);
    const double water_vol = 3.0;
    FillDepressions(pit_cell, dep_labels, water_vol, topo, label, wtd, ocean_level);

    const auto W = 1.5;
    const Array2D<double> wtd_good = {
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
        { 0,  0,  0,  0,  W,  W,  0,  0,  0,  0},
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
    };

    CHECK(MaxArrayDiff(wtd,wtd_good)<1e-6);
  }

  SUBCASE("Sill Cell Aborbs some water"){
    wtd(4,3) = -1;
    const double water_vol = 5.0;
    FillDepressions(pit_cell, dep_labels, water_vol, topo, label, wtd, ocean_level);

    const auto W = 2.0;
    const Array2D<double> wtd_good = {
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
        { 0,  0,  0,  0,  W,  W,  0,  0,  0,  0},
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
    };

    CHECK(MaxArrayDiff(wtd,wtd_good)<1e-6);
  }

  SUBCASE("Passes over a saddle"){
    const double water_vol = 19.0;
    FillDepressions(pit_cell, dep_labels, water_vol, topo, label, wtd, ocean_level);

    const Array2D<double> wtd_good = {
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
        { 0,  0,  0,  0,  3,  3,  0,  0,  0,  0},
        { 0,  0,  0,  0,  1,  1,  0,  0,  0,  0},
        { 0,  0,  0,  0,  2,  3,  0,  0,  0,  0},
        { 0,  0,  0,  0,  3,  3,  0,  0,  0,  0},
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
    };

    CHECK(MaxArrayDiff(wtd,wtd_good)<1e-6);
  }

  SUBCASE("Passes over a saddle and sill absorbs some"){
    const double water_vol = 19.5;
    wtd(6,4) = -1;
    FillDepressions(pit_cell, dep_labels, water_vol, topo, label, wtd, ocean_level);

    const Array2D<double> wtd_good = {
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
        { 0,  0,  0,  0,  3,  3,  0,  0,  0,  0},
        { 0,  0,  0,  0,  1,  1,  0,  0,  0,  0},
        { 0,  0,  0,  0,  2,  3,-.5,  0,  0,  0},
        { 0,  0,  0,  0,  3,  3,  0,  0,  0,  0},
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
    };

    CHECK(MaxArrayDiff(wtd,wtd_good)<1e-6);
  }
}



TEST_CASE("Randomized Heavy Flooding vs Priority-Flood"){
  std::mt19937_64 gen;
  std::uniform_int_distribution<int>      size_dist(30,500);
  std::uniform_int_distribution<uint32_t> seed_dist;

  #pragma omp parallel for
  for(int i=0;i<500;i++){
    int      this_size;
    uint32_t this_seed;
    #pragma omp critical
    {
      this_size = size_dist(gen);
      this_seed = seed_dist(gen);
    }
    auto dem = perlin(this_size, this_seed);

    Array2D<dh_label_t> labels  (dem.width(), dem.height(), NO_DEP );
    Array2D<flowdir_t>  flowdirs(dem.width(), dem.height(), NO_FLOW);

    for(int y=0;y<dem.height();y++)
    for(int x=0;x<dem.width(); x++){
      if(dem.isEdgeCell(x,y)){
        dem(x,y)    = -1;
        labels(x,y) = OCEAN;
      }
    }

    auto deps = GetDepressionHierarchy<double,Topology::D8>(dem, labels, flowdirs);

    //wtd with a *lot* of initial surface water
    Array2D<double> wtd(dem.width(), dem.height(), 100);

    FillSpillMerge(dem, labels, flowdirs, deps, wtd, -1.0);

    for(auto i=dem.i0(); i<dem.size(); i++){
      if(!dem.isNoData(i))
        dem(i) += wtd(i);
    }

    auto comparison_dem = dem;
    PriorityFlood_Zhou2016(comparison_dem);

    CHECK(MaxArrayDiff(comparison_dem,dem)<1e-6);
  }
}
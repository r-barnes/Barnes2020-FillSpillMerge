#include "doctest.h"
#include <fsm/fill_spill_merge.hpp>

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

  CHECK(wtd==wtd_good);
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